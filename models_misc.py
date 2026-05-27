import numpy as np

__all__ = [
    'benyahia2005', 'ribeiro2010', 'gas_props', 'water_props', 'liquid_holup', 'calculate_G_cond', 'massHTU', 'outputGasComposition'
  ]

############################################################################################
def benyahia2005(D, D_p_eff, phi_s, **kwargs): # void fraction
  """
  Calculate mean void fraction of a packed bed using Benyahia & O'Neill equation
  https://doi.org/10.1080/02726350590922242
  """
  D_p = D_p_eff
  ε = ε_eff
  valid = (
    (D/D_p > 1.5) & (D/D_p < 50)            # benyahia2005, p. 10
    & (phi_s > 0.42)                        # benyahia2005, p. 10
    
  )

  ε = 0.1504 + 0.2024 / phi_s + 1.0814 / (D / D_p + 0.1226)**2  # mean void fraction, eq. (4)
  ε_array = np.where(valid, ε, np.nan)

  return { 'ε': ε_array }

############################################################################################
def ribeiro2010(D, D_p_eff, phi_s, **kwargs): # void fraction
  """
  Calculate mean void fraction of a packed bed using Ribeiro et al. equation
  https://doi.org/10.1080/02726350590922242
  """
  D_p = D_p_eff
  ε = ε_eff
  valid = (
    (phi_s > 0.95) &                        # assumed
    (D/D_p >= 2) & (D/D_p <= 19)            # ribeiro2010
  )

  ε = 0.373 + 0.917 + np.exp(-0.824 * D / D_p)  # mean void fraction
  ε_array = np.where(valid, ε, np.nan)

  return { 'ε': ε_array }

############################################################################################
def gas_props(M_CO2, M_H2O, M_O2, T, P=101325, props=['ρ', 'μ', 'κ', 'Cp', 'MW']):
  """
  Vectorized property calculator for flue gas mixtures.

  Args:
    M_...: np.arrays of mole fractions for components.
    T (np.array): Temperature in Degrees Celcius.
    P (float/np.array): Pressure in Pascals.
    props (list): Properties to return ['ρ', 'μ', 'Cp'].

  Returns:
    dict: Requested properties as NumPy arrays.
  """
  # Defensive Conversion to NumPy Arrays
  M_CO2 = np.asarray(M_CO2)
  M_H2O = np.asarray(M_H2O)
  M_O2 = np.asarray(M_O2)
  T = np.asarray(T)
  P = np.asarray(P) # Handles both scalar 101325 and arrays

  R = 8.314
  T_K = T + 273.15
  t = T_K / 1000.0
  M_N2 = 1.0 - M_CO2 - M_H2O - M_O2
  components = {'CO2': M_CO2, 'H2O': M_H2O, 'O2': M_O2, 'N2': M_N2}

  # Shomate Coefficients [A, B, C, D, E] for J/mol*K
  shomate = {
    'N2':  [26.092, 8.2188, -1.9761, 0.1592, 0.0444],
    'H2O': [30.092, 6.8325, 6.7934, -2.5345, 0.0821],
    'CO2': [24.997, 55.186, -33.691, 7.9483, -0.1366],
    'O2':  [29.659, 6.1372, -1.1865, 0.0957, -0.2196]
  }

  # Species data: Molar Mass (kg/mol), Sutherland C, mu0 (Pa-s), T0 (K), k0 (W/m-K)
  species_data = {
    'N2':  {'MW': 0.02801, 'C': 111, 'mu0': 1.781e-5, 'T0': 288, 'k0': 0.024},
    'H2O': {'MW': 0.01802, 'C': 961, 'mu0': 1.12e-5,  'T0': 350, 'k0': 0.018},
    'CO2': {'MW': 0.04401, 'C': 240, 'mu0': 1.48e-5,  'T0': 293, 'k0': 0.014},
    'O2':  {'MW': 0.03199, 'C': 127, 'mu0': 2.018e-5, 'T0': 292, 'k0': 0.024}
  }

  results = {}
  mix_MW = sum(components[gas] * species_data[gas]['MW'] for gas in components)
  results['MW'] = mix_MW

  # 1. Density (Ideal Gas Law)
  if 'ρ' in props:
    results['ρ'] = (P * mix_MW) / (R * T_K)

  # 2. Viscosity (Wilke Mixing Rule)
  if 'μ' in props:
    # Calculate pure component viscosities using Sutherland's Law
    mus = {g: species_data[g]['mu0'] * ((species_data[g]['T0'] + species_data[g]['C']) /
            (T_K + species_data[g]['C'])) * (T_K / species_data[g]['T0'])**1.5
            for g in components}

    mix_mu = np.zeros_like(T_K, dtype=float)
    for i in components:
      phi_sum = 0
      for j in components:
        # Wilke interaction parameter
        num = (1 + (mus[i]/mus[j])**0.5 * (species_data[j]['MW']/species_data[i]['MW'])**0.25)**2
        den = (8 * (1 + species_data[i]['MW']/species_data[j]['MW']))**0.5
        phi_sum += components[j] * (num / den)
      mix_mu += (components[i] * mus[i]) / phi_sum
    results['μ'] = mix_mu

  # 3. Thermal Conductivity (Simplified T-scaling)
  if 'κ' in props:
    results['κ'] = sum(components[g] * species_data[g]['k0'] * (T_K / 273.15)**1.03 for g in components)

  # 4. Specific Heat (Cp) - Mass Basis
  mix_cp_molar = np.zeros_like(T_K)
  for g, frac in components.items():
    A, B, C, D, E = shomate[g]
    cp_species = A + B*t + C*(t**2) + D*(t**3) + E/(t**2)
    mix_cp_molar += frac * cp_species
  results['Cp'] = mix_cp_molar / mix_MW  # Convert J/mol-K to J/kg-K

  return [results[p] for p in props]

############################################################################################
def water_props(T=20, P=101325, props=['ρ', 'μ', 'κ', 'Cp']):
  """
  Vectorized properties for liquid water at 1 atm (101325 Pa).
  Accuracy range: 10 - 100 °C.
  """
  T = np.asarray(T, dtype=float)
  P = np.asarray(P, dtype=float)
  results = {}

  if 'ρ' in props:
    # Density (Kell correlation, kg/m3)
    a = [999.83311, 0.0752, -0.0089, 7.36413e-5, -4.74639e-7, 1.34888e-9]
    results['ρ'] = (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4 + a[5]*T**5)

  if 'μ' in props:
    # Viscosity (Vogel-style correlation, Pa-s)
    T_K = T + 273.15
    results['μ'] = 0.001 * np.exp(-3.7188 + 578.919 / (T_K - 137.546))

  if 'κ' in props:
    # Thermal Conductivity (W/m-K)
    results['κ'] = 0.561 + 0.0019 * T - 8.0e-6 * T**2

  if 'Cp' in props:
    # Specific Heat Cp (J/kg-K)
    results['Cp'] = 4217.4 - 3.72 * T + 0.141 * T**2 - 0.00227 * T**3 + 1.28e-5 * T**4

  return [results[p] for p in props]

############################################################################################
def liquid_holup(Ratio_L__G, D, D_p, phi_s, ε, U_s_G_in, ρ_G_avg, ρ_L_avg, μ_L_avg, G_cond, **kwargs):
  """
  Quick estimate of maximum possible condensate (perfect column, infinitely tall). No packing data needed.
  """
  # Defensive Conversion to NumPy Arrays
  Ratio_L__G = np.asarray(Ratio_L__G)
  D          = np.asarray(D)
  D_p        = np.asarray(D_p)
  phi_s      = np.asarray(phi_s)
  ε          = np.asarray(ε)
  U_s_G_in   = np.asarray(U_s_G_in)
  ρ_G_avg    = np.asarray(ρ_G_avg)
  ρ_L_avg    = np.asarray(ρ_L_avg)
  μ_L_avg    = np.asarray(μ_L_avg)
  G_cond     = np.asarray(G_cond)

  G_fg_in = U_s_G_in * ρ_G_avg * np.pi * D**2 / 4
  G_wtr_in = G_fg_in * Ratio_L__G
  G_L = G_wtr_in + G_cond

  U_s_L = G_L / ρ_L_avg / np.pi / D**2 * 4 # Superficial liquid (rinse water + condensate) velocity
  A_ss = 6 * (1 - ε) / phi_s / D_p # specific surface area [m²/m³]
  Re_L = (ρ_L_avg * U_s_L * D_p) / (μ_L_avg * (1 - ε)) # Reynolds Number (Liquid)

  # Inertial regime (Re_L < Re_L_cr) https://doi.org/10.1016/j.cjche.2023.06.009 https://doi.org/10.1016/0950-4214(89)80016-7
  Fr_L = U_s_L**2 * A_ss / 9.81 / ε**4.65 # liquid Froude number
  h_L_inertial = 0.555 * Fr_L**0.333

  # Viscous regime (Re_L > Re_L_cr): Billet and Schultes equation https://doi.org/10.1205/026387699526520
  h_L_viscous = (12 / 9.81 * μ_L_avg / ρ_L_avg * U_s_L * A_ss**2)**0.333

  # Sigmoid Blend
  k = 0.2 # blending factor
  Re_L_cr = 60 # critical Reynolds
  exponent = np.clip(k * (Re_L - Re_L_cr), -100, 100) # Clip the exponent to prevent overflow
  w_v = 1 / (1 + np.exp(exponent))
  h_L = (w_v * h_L_viscous) + ((1 - w_v) * h_L_inertial)

  return h_L # the fraction of total bed volume occupied by liquid, dimensionless

############################################################################################
def calculate_G_cond(M_G_in_H2O, M_G_in_CO2, M_G_in_O2, U_s_G_in, T_G_in, T_wtr_in, D, P=101325):
  # Defensive Conversion to NumPy Arrays
  M_G_in_H2O = np.asarray(M_G_in_H2O)
  M_G_in_CO2 = np.asarray(M_G_in_CO2)
  M_G_in_O2  = np.asarray(M_G_in_O2)
  U_s_G_in   = np.asarray(U_s_G_in)
  T_G_in     = np.asarray(T_G_in)
  T_wtr_in   = np.asarray(T_wtr_in)
  D          = np.asarray(D)
  P          = np.asarray(P)

  # 1. Define Species MWs
  MW = {'H2O': 0.01802, 'CO2': 0.04401, 'O2': 0.03199, 'N2': 0.02801}
  M_G_in_N2 = 1.0 - M_G_in_H2O - M_G_in_CO2 - M_G_in_O2

  # 2. Get Inlet Dry Gas Mass Flow
  ρ_G_in, = gas_props(M_G_in_H2O, M_G_in_CO2, M_G_in_O2, T_G_in, props=['ρ'])
  G_total_in = U_s_G_in * ρ_G_in * (np.pi * D**2 / 4)

  # Mass fraction of water in inlet
  Y_H2O_in = (M_G_in_H2O * MW['H2O']) / (
      M_G_in_H2O*MW['H2O'] + M_G_in_CO2*MW['CO2'] + M_G_in_O2*MW['O2'] + M_G_in_N2*MW['N2']
  )
  m_dry = G_total_in * (1 - Y_H2O_in) # Constant dry gas mass flow

  # 3. Calculate Inlet Humidity Ratio (omega_in)
  MW_dry_in = (M_G_in_CO2*MW['CO2'] + M_G_in_O2*MW['O2'] + M_G_in_N2*MW['N2']) / (1 - M_G_in_H2O)
  omega_in = (M_G_in_H2O * MW['H2O']) / ((1 - M_G_in_H2O) * MW_dry_in)

  # 4. Estimate Outlet Humidity Ratio (omega_out)
  # Assumption: Gas exits saturated (RH=100%) at water inlet temperature
  # P_sat from Antoine or similar
  P_sat_wtr = 10** (8.07131 - 1730.63 / (T_wtr_in + 233.426)) * 133.322 # Pa
  M_H2O_out = P_sat_wtr / P

  # Note: Dry gas composition ratio (CO2/N2/O2) is assumed constant
  omega_out = (M_H2O_out * MW['H2O']) / ((1 - M_H2O_out) * MW_dry_in)

  # 5. Condensate Rate
  G_cond = m_dry * (omega_in - omega_out)
  return np.maximum(0, G_cond) # Condensation only if omega_in > omega_out

############################################################################################
def massHTU(flux, k, α_e, C):
  """
  Calculate mass Height of Transfer Unit.

  flux is a molar flux [kmol/(m²·s)]
  C is a molar concentration [kmol/m³]
  k is a mass transfer coeff [m/s]
  α_e is an effective interfacial area [1/m]
  """
  # Defensive Conversion to NumPy Arrays
  flux = np.asarray(flux)
  C   = np.asarray(C)
  k   = np.asarray(k)
  α_e = np.asarray(α_e)

  if k is not None and α_e is not None:
    return flux / (k * α_e * C)
  else:
    return None

############################################################################################
def outputGasComposition(M_G_in_CO2, M_G_in_H2O, M_G_in_O2, Ratio_L__G, U_s_G_in, T_G_in, T_wtr_in, D, Z, P=101325, k_G=None, k_L=None, α_e=None, HTU_G=None, HTU_L=None):
  """
  Calculate outlet gas composition after condensation using molar fluxes.
  Realistic prediction for a given packed height and mass transfer coefficients. Requires packing parameters.
  """
  # Convert to arrays
  M_G_in_H2O = np.asarray(M_G_in_H2O)
  M_G_in_CO2 = np.asarray(M_G_in_CO2)
  M_G_in_O2  = np.asarray(M_G_in_O2)
  Ratio_L__G  = np.asarray(Ratio_L__G)
  U_s_G_in = np.asarray(U_s_G_in)
  T_G_in   = np.asarray(T_G_in)
  T_wtr_in = np.asarray(T_wtr_in)
  D        = np.asarray(D)
  Z        = np.asarray(Z)
  P        = np.asarray(P)

  A = np.pi * D**2 / 4.0                     # cross-sectional area (m²)

  # ---------- 1. Inlet gas mass and molar flows ----------
  # Get gas density (kg/m³) and average molecular weight (kg/kmol)
  rho_G_in, MW_G_in = gas_props(M_H2O=M_G_in_H2O, M_CO2=M_G_in_CO2, M_O2=M_G_in_O2, T=T_G_in, P=P, props=['ρ', 'MW'])
  gas_mass_flow_in = U_s_G_in * rho_G_in * A          # kg/s
  G_total_in_kmol_s = gas_mass_flow_in / MW_G_in      # kmol/s

  # Inert (dry gas) molar flow – constant
  G_inert_kmol_s = G_total_in_kmol_s * (1 - M_G_in_H2O)   # kmol/s

  # ---------- 2. Inlet liquid molar flow from mass ratio ----------
  liquid_mass_flow_in = R_LG * gas_mass_flow_in         # kg/s
  L_in_kmol_s = liquid_mass_flow_in / 0.01802           # kmol/s

  # ---------- 3. Equilibrium at water temperature ----------
  P_sat_mmHg = 10**(8.07131 - 1730.63 / (T_wtr_in + 233.426))
  P_sat = P_sat_mmHg * 133.322            # Pa
  y_star = P_sat / P                      # equilibrium mole fraction
  m = y_star                              # slope of equilibrium line (x=1)

  C_L = 55.5                               # liquid molar density (kmol/m³)
  R_gas = 8314.0                           # J/(kmol·K)

  # ---------- 4. Initial guess for outlet water mole fraction ----------
  y_out = y_star.copy()
  max_iter = 10

  for iteration in range(max_iter):
    # ----- 4a. Material balances -----
    G_total_out_kmol_s = G_inert_kmol_s / (1 - y_out)
    delta_water = G_total_in_kmol_s * M_G_in_H2O - G_total_out_kmol_s * y_out
    L_out_kmol_s = L_in_kmol_s + delta_water

    # Average molar fluxes (kmol·m⁻²·s⁻¹)
    G_avg = (G_total_in_kmol_s + G_total_out_kmol_s) / (2 * A)
    L_avg = (L_in_kmol_s + L_out_kmol_s) / (2 * A)

    # ----- 4b. Average gas properties -----
    T_avg = (T_G_in + T_wtr_in)/2
    C_total_avg = P / (R_gas * (T_avg + 273.15))   # total gas molar concentration (kmol/m³)

    # ----- 4c. HTU calculations (molar flux basis) -----
    HTU_G = massHTU(flux=G_avg, k=k_G, α_e=α_e, C=C_total_avg) if HTU_G is None else HTU_G
    HTU_L = massHTU(flux=L_avg, k=k_L, α_e=α_e, C=55.5) if HTU_L is None else HTU_L
    if HTU_G is None or HTU_L is None:
      return None

    stripping = m * G_avg / L_avg
    HTU_OG = HTU_G + stripping * HTU_L

    # ----- 4d. NTU and outlet mole fraction -----
    NTU_OG = Z / HTU_OG
    y_out_new = y_star + (M_G_in_H2O - y_star) * np.exp(-NTU_OG)

    if np.all(np.abs(y_out_new - y_out) < 1e-8):
      y_out = y_out_new
      break
    y_out = y_out_new

  # ---------- 5. Final outputs ----------
  G_total_out_kmol_s = G_inert_kmol_s / (1 - y_out)
  G_L_out = (G_total_in_kmol_s * M_G_in_H2O - G_total_out_kmol_s * y_out) * 18.02

  dry_inert_sum = 1 - M_G_in_H2O
  M_G_out_CO2 = M_G_in_CO2 * dry_inert_sum / (1 - y_out)
  M_G_out_O2  = M_G_in_O2  * dry_inert_sum / (1 - y_out)
  MW_G_out = gas_props(M_H2O=y_out, M_CO2=M_G_out_CO2, M_O2=M_G_out_O2, P=P, props=['MW'])

  return {
      'M_G_out_H2O': y_out,
      'G_G_out': G_total_out_kmol_s * MW_G_out,
      'G_L_out': G_L_out, # condensate output in kg/s
      'M_G_out_CO2': M_G_out_CO2,
      'M_G_out_O2': M_G_out_O2,
  }