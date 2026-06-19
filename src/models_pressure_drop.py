import numpy as np

__all__ = [
    'carman1937', 'morcom1946', 'rose1949', 'leva1949', 'ergun1952',
    'wentz1963', 'kurten1966', 'handley1968',
    'mehta1969', 'hicks1970', 'tallmadge1970', 'reichelt1972', 'kuo1978',
    'macdonald1979', 'kta1981', 'jones1983', 'foscolo1983', 'fahien1983', 'meyer1985',
    'paterson1986', 'stichlmair1989', 'watanabe1989', 'comiti1989', 'fand1990',
    'foumeny1993', 'lee1994', 'liu1994', 'hayes1995', 'avontuur1996', 'critoph1996',
    'oneill1997', 'raichura1999', 'eisfeld2001', 'yu2002', 'felice2004',
    'nemec2005', 'montillet2007', 'carpinlioglu2008', 'ozahi2008', 'cheng2011',
    'harrison2013', 'erdim2015', 'guo2017', 'seckendorff2020', 'cheng2021', 'reger2023',
    'dixon2023', 'wu2025', 'xu2026',
]

############################################################################################
def test(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using test equation
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
     (phi_s >= 0.98)        # test
  )

  F_s = 150 / Re_m + 1.75
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def carman1937(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Carman--Kozeny equation
  https://doi.org/10.1016/S0263-8762(97)80003-2
  """

  D_p = D_p_eff
  ε = ε_eff
  S = 6 * (1 - ε) / phi_s / D_p           # specific surface area of particles [m²/m³]
  S1 = S + 4 / D                          # wall correction
  Re_1 = Re_p * phi_s / 6 / (1 - ε)
  valid = (
    (Re_1 >= 0.01) & (Re_1 <= 10000)   # carman1937, Summary p.15
    & (D/D_p >= 2)                     # carman1937, deduced from Appendix III and wall-effect discussion
    & (phi_s >= 0.95)                  # carman1937, correlation for spherical particles (Eq. 9a, Summary)
    & (ε >= 0.30) & (ε <= 0.90)        # carman1937, from Fig.1 and Table IX
  )

  F_s = 1 / Re_1+ 0.4 / Re_1**0.1
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 * S1 / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def morcom1946(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Morcom equation
  https://scholar.google.com/scholar_lookup?title=Fluid%20flow%20through%20granular%20materials&publication_year=1946&author=A.R.%20Morcom
  A.R. Morcom, Fluid flow through granular materials, Chemical Engineering Research and Design 24 (1946), pp. 30-43.
  no text available, reference only: erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  """

  D_p = D_p_eff
  ε = ε_eff
  valid = (
    (Re_p < 750)                         # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
    & (phi_s >= 0.95)                    # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
    & (D/D_p > 5)                        # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  )

  F_p = (784.8 / Re_p + 13.73) * (0.405 / ε)**3
  ΔP = L * F_p * ρ_G_avg * U_s_G_in**2 / D_p  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def rose1949(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Rose & Rizk equation
  https://doi.org/10.1243/PIME_PROC_1949_160_047_02
  """

  D_p = D_p_eff
  ε = ε_eff
  valid = (
    (Re_p >= 100) & (Re_p <= 100000)     # rose1949, Fig. 16
    & (D/D_p >= 1)                       # rose1949, Fig. 5
    & (phi_s >= 0.3)                     # rose1949, deduced from Table 1
    & (ε >= 0.3) & (ε <= 0.9)            # rose1949 from Abstract
  )

  f = 36 * np.exp(-0.0915 * ε) + 0.055  # from Fig. 1
  F_p = (1000 / Re_p + 125 / Re_p**0.5 + 14) * f
  ΔP = L * F_p * ρ_G_avg * U_s_G_in**2 / D_p  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def leva1949(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Leva equation
  D.W. Green, R.H. Perry (Eds.), Perry's Chemical Engineers' Handbook (Eight edition), Chem. Eng., 56 (1949), pp. 115-117

  !! didn't find orig paper, used approximations proposed by erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  """

  D_p = D_p_eff
  ε = ε_eff
  valid = (
    (Re_p < 10000)                       # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
    & (phi_s >= 0.95)                    # assumed
    #& (ε >= 0.36) & (ε <= 0.4)           # assumed
  )

  x = np.log10(Re_p)
  S1 = 7.60657 - 19.2986 * x + 21.02695 * x**2  - 10.96663 * x**3 + 3.02928 * x**4 - 0.42867 * x**5 + 0.02453 * x**6
  n = np.where(Re_p < 11.5, 1, S1)
  S2 = 1.982535 - 1.0218594 * Re_p + 0.0295464 * Re_p**2 + 0.0269893 * Re_p**3 + 0.0024996 * Re_p**4 - 0.0008754 * Re_p**5
  F_m = 10**S2
  ΔP = 2 * L * F_m * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε)**(3 - n) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def ergun1952(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Ergun equation
  S. Ergun, Fluid flow through packed columns, Chem. Eng. Progress 48(2) (1952), pp. 89-94.
  https://api.semanticscholar.org/CorpusID:135806248
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m >= 1.2) & (Re_m <= 4200)   # erdim2015revisit https://doi.org/10.1016/j.ijheatmasstransfer.2024.126620, range Re₁ = 0.2-700, where Re₁ = Re_m/6
    & (D/D_p > 10)                   # ergun1952, from page 5
    & (phi_s >= 0.95)                # ergun1952
  )

  F_s = 150 / Re_m + 1.75
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def brauer1957removed(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Brauer equation

  !!! Brauer doesn't state F_s = 160 / Re_m + 3.1 / Re_m**0.1 in his 1957/1960 papers, but this correlation is quite popular among reviwers

  1957 https://doi.org/10.1002/cite.330291207 — single-phase flow
  1960 https://doi.org/10.1002/cite.330320905 — two-phase flow

  his models are for hollow Rachig rings only. It's safer to skip.
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m >= 10) & (Re_m <= 20000)    # brauer1957, Fig. 5
     & (L/D >= 0.86) & (L/D <= 3.59)  # brauer1957, from page 5: H = 30-80 cm, D = 22.3 & 35 cm
     & (D/D_p >= 8) & (D/D_p <= 63)   # brauer1957, from Table 1 and column diameters
     & (phi_s > 0.85) & (phi_s < 0.89)# brauer1957, deduced for equilateral cylinders (Raschig rings)
  )

  ## F_s = 160 / Re_m + 3.1 / Re_m**0.1  Brauer doesn't state this
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def wentz1963(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Wentz & Thodos equation
  https://doi.org/10.1002/aic.690090118
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m > 2550) & (Re_m < 64900)        # wentz1963, from abstract
    & (L/D_p >= 5)                        # wentz1963, "five layers of spheres"
    & (D/D_p >= 11)                       # wentz1963, verified for D/D_p = 11.38 from D=14 in, D_p=1.23 in (page 2)
    & (phi_s >= 0.95)                     # wentz1963
    & (ε > 0.354) & (ε < 0.882)           # wentz1963, from Table I and void fraction description (page 1-2)
    
  )
  A = 0.396 - 0.045 * np.exp(-0.47 * (L / D_p - 5)) # my assumption of asymptotic approximation (A=0.351 at 5 layers and growing to 0.396)
  F_s = A / (Re_m**0.05 - 1.2)
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def kurten1966(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Kürten et al. equation
  https://doi.org/10.1002/cite.330380905

  !! orig article states c_w only and for a single partcle only. extrapolation to packed bed made by Ozahi.
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p >= 0.1) & (Re_p <= 4000)       # kurten1966
    & (phi_s >= 0.95)                    # kurten1966
  )
  c_w = (21 / Re_p + 6 / Re_p**0.5 + 0.28) # kurten1966
  F_s = 25 / 4 * (1 - ε) * c_w             # ozahi2008 https://doi.org/10.1163/156855208X314985
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def handley1968(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Handley & Heggs equation
  https://scholar.google.com/scholar?cluster=8679261533528586850&hl=en&as_sdt=0,5
  D. Handly, Momentum and heat transfer mechanisms in regular shaped packings. Trans. Inst. Chem. Eng., 46 (1968), pp. 251-259.
  !! no orig paper found
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p > 200) & (Re_p < 13000)    # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
    & (phi_s >= 0.8)                 # cylinders are mentioned by critoph1996, p. 2
    #& (ε >= 0.38) & (ε <= 0.42)      # assumed
  )

  F_s = 368 / Re_m + 1.24
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def mehta1969(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Mehta & Hawley equation
  https://doi.org/10.1021/i260030a021
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  M = 1 + 0.66667 * D_p / (1 - ε) / D  # Eq. (4)
  valid = (
    (Re_m/M <= 10)                   # mehta1969, Fig. 3
    & (D/D_p >= 7)                   # mehta1969 , from Experimental
    & (phi_s >= 0.95)                # mehta1969, spherical glass beads, p.1
  )

  F_s = 150 * M**2 / Re_m + 1.75 * M
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def hicks1970(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Hicks equation
  https://doi.org/10.1021/i160035a032
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m > 300) & (Re_m < 60000)    # hicks1970, p. 3
    & (phi_s >= 0.95)                # hicks1970
  )

  F_s = 6.8 * Re_m**1.2
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def tallmadge1970(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Tallmadge equation
  https://doi.org/10.1002/aic.690160639
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m > 0.1) & (Re_m < 100_000)  # tallmadge1970, Eq. (2)
    & (phi_s >= 0.95)                # tallmadge1970
    & (ε >= 0.35) & (ε <= 0.88)      # tallmadge1970, porosity range 35–88% from Wentz & Thodos data
  )

  F_s = 150 / Re_m + 4.2 / Re_m**0.166667
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def reichelt1972(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Reichelt equation
  https://doi.org/10.1002/cite.330441806
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  Re_w = Re_m / (1 + 0.67 * D_p / D / (1 - ε))
  valid = (
    (Re_w >= 0.2) & (Re_w <= 30000)      # reichelt1972
    & (D/D_p > 1.7)                      # reichelt1972
    & (phi_s >= 0.85)                    # reichelt1972 implicitly: around 0.87 for equilateral cylinders and >0.95 for spheres
  )

  B_sphere = D**4 / D_p**4 / (1.5 + 0.88 * D**2/D_p**2)**2 # 1/√K₂ = 1.5 / (D/D_p)² + 0.88, eq. (15)
  B_cylinder = D**4 / D_p**4 / (2 + 0.8 * D**2/D_p**2)**2 # 1/√K₂ = 2 / (D/D_p)² + 0.8, eq. (17)
  F_s = np.where(phi_s > 0.9, 150 / Re_w + B_sphere, 200 / Re_w + B_cylinder)
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def kuo1978(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Kuo & Nydegger equation
  K.K. Kuo, C.C. Nydegger, Flow resistance measurements and correlation in a packed bed of WC 870 ball propellants. J. of Ballistics 2(1) (1978), pp. 1-25.

  !! no orig paper found, but have confirmation in two sources ( erdim2015revisit and report to Defense Technical Information Center (DTIC, https://apps.dtic.mil/sti/tr/pdf/ADA091300.pdf))
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p >= 460) & (Re_p <= 14600)        # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
    #& (phi_s >= 0.8) & (phi_s < 0.95)      # assumed as the exp particles were not spherical
    & (ε >= 0.376) & (ε <= 0.39)           # DTIC report
  )

  F_s = 276.23 / Re_m + 5.05 / Re_m**0.13          # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def macdonald1979(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Macdonald et al. equation
  https://doi.org/10.1021/i160071a001
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m <= 10000)                        # macdonald1979, Fig. 10
    & (phi_s >=0.6) & (phi_s <= 1)         # macdonald1979, spherical φ_s=1 (Gupte data), irregular φ_s=0.6 (p.6, Table III), equilateral cylinders ~0.87
    & (ε >= 0.36) & (ε <= 0.92)            # macdonald1979, Conclusion 9
  )

  F_s = 180 / Re_m + 1.8
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def kta1981(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using KTA equation
  https://www.kta-gs.de/e/standards/3100/3102_3_engl_1981_03.pdf
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  cond1 = (Re_m <= 1000) & (D/D_p > 34.35 - 0.06 * Re_m)
  cond2 = (Re_m > 1000) & (D/D_p > 5)
  valid = (
    (cond1 | cond2)                 # kta1981
    & (phi_s >= 0.95)               # kta1981
    & (ε >= 0.36) & (ε <= 0.42)     # kta1981
    & (L/D_p > 5)                   # kta1981
  )

  F_s = 160 / Re_m + 3 / Re_m**0.1
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def jones1983(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Jones & Krier equation
  https://doi.org/10.1115/1.3240959
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m > 733) & (Re_m < 126670)             # jones1983, Table 2
    & (L/D > 3.5)                              # jones1983, verified for L/D=3.94 from test section (page 2)
    & (D/D_p >= 20)                            # jones1983, min from Dc/Db>=20, (page 3, Table A1)
    & (phi_s >= 0.95)                          # jones1983, spherical glass beads (page 1)
    & (ε >= 0.372) & (ε <= 0.436)              # jones1983, Table A1 and Table 1
  )

  F_s = 150 / Re_m + 3.89 / Re_m**0.13
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def foscolo1983(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Foscolo et al. equation
  https://doi.org/10.1016/0009-2509(83)80045-1
  """

  D_p = D_p_eff
  ε = ε_eff
  valid = (
    (phi_s >= 0.95)              # foscolo1983, spherical particles assumed in Intoroduction
    & (ε >= 0.4)                 # foscolo1983, Figs. 2-6
  )

  F_p = (17.3 / Re_p + 0.336) * (1 - ε) / ε**4.8
  ΔP = L * F_p * ρ_G_avg * U_s_G_in**2 / D_p  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def fahien1983(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Fahien & Schriver equation
  Fahien, R. W., & Schriver, C. B. (1983). Paper presented and Denver meeting of AIChE. Fundamentals of Transport Phenomena, McGraw-Hill, New York.
  !! no original text found
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (phi_s >= 0.95)                    # guessed
    #& (ε >= 0.36) & (ε <= 0.42)        # guessed
    #& (D/D_p > 5)                      # guessed
  )

  q = np.exp(-ε**2 * (1 - ε) * Re_m / 12.6)
  f_1L = 136 / (1- ε)**0.38
  f_1T = 29 / (1 - ε)**1.45 / ε**2
  f2 = 1.87 * ε**0.75 / (1 - ε)**0.26
  F_s = q * f_1L / Re_m + (1 - q) * (f2 + f_1T / Re_m) # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def meyer1985(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Meyer & Smith equation
  https://doi.org/10.1021/i100019a013

  Eq. (22) accounts for particle roughness q' via channel roughness q but the coeffs 0.49 and 0.8 
  were derived for a special case when q/D_p_eff ≈ 0.15; thus, we use the model without roughness.
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  S_v = 6 * (1 - ε) / phi_s / D_p    # specific surface area of particles [m²/m³]
  G0 = U_s_G_in * ρ_G_avg            # Superficial mass flux [kg/(m²·s)]
  Re_s = Re_m * phi_s / 6            # Re_s = G0 / μ / S_v = ρ·U_s·D_p·ϕ_s / (6μ·(1-ε))
  valid = (
    (Re_s >= 0.01) & (Re_s <= 1000)  # meyer1985, Summary
    & (ε >= 0.18) & (ε <= 0.67)      # meyer1985, Abstract
  )

  
  ΔP = L * G0**2 * S_v / ρ_G_avg / ε**4.1 * (2.5 * 6 / Re_m / phi_s + 0.077)  # Total pressure drop, Eq. (22) without q
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def paterson1986(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Paterson et al. equation
  Paterson, W. R., Burns, J. R. M., Griffiths, N. B., Kesterton, K. R., & Paveley, A. J. (1986). Experimental studies of transport processes in packed beds of low tube-to-particle diameter ratio. In World Congress III of Chemical Engineering Tokyo (pp. 304-307).

  !! no orig paper found
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p > 25) & (Re_p < 900)           # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
    & (phi_s >= 0.95)                    # assumed
    #& (ε >= 0.3) & (ε <= 0.4)            # assumed
    & (D/D_p > 3.5) & (D/D_p < 22)       # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  )

  A = 1 + 1.22 * D_p / D
  B = np.exp(1.66 * ((1 - D_p/D)**2 - 1))
  F_s = 150 * A / Re_m + 1.75 * B     # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def stichlmair1989(D, L, D_p_eff, ε_eff, h_L, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Stichlmair equation
  https://doi.org/10.1016/0950-4214(89)80016-7
  """

  D_p = D_p_eff
  ε = ε_eff
  μ_G_avg = U_s_G_in * D_p * ρ_G_avg / Re_p
  valid = (
    (Re_p >= 0.01) & (Re_p <= 100_000)      # stichlmair1989, Figs. 1, 2
    & (phi_s >= 0.95)                       # stichlmair1989
    & (μ_G_avg <= 0.006)                    # stichlmair1989, Fig. 7
  )
  F_s = 24 / Re_p + 4 / Re_p**0.5 + 0.4       # for spherical particles only (Fig. 2)
  c = -24 / Re_p / F_s - 2 / Re_p**0.5 / F_s  # for spherical particles only (Fig. 2)
  ΔP_dry = 0.75 * L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**4.65  # Total pressure drop (dry)
  ΔP = ΔP_dry * ((1 - ε * (1 - h_L / ε)) / (1 - ε))**(2/3 + c/3) * (1 - h_L/ε)**-4.65
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def watanabe1989(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Watanabe equation
  H. Watanabe. Drag coefficient and voidage function on fluid-flow through granular packed-beds. Int. J. of Eng. Fluid Mechanics 2(1) (1989), pp. 93-108.

  !! no orig paper found
  """

  D_p = D_p_eff
  ε = ε_eff
  valid = (
    (Re_p > 0.1) & (Re_p < 4000)         # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
    & (phi_s >= 0.95)                     # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
    #& (ε >= 0.3) & (ε <= 0.4)            # assumed
    #& (D/D_p > 10)                       # assumed
  )

  F_s = 6.25 * (21 / Re_p + 6 / Re_p**0.5 + 0.28) * (1 - ε)     # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def comiti1989(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Comiti & Renaud equation
  https://doi.org/10.1016/0009-2509(89)80031-4
  """

  D_p = D_p_eff
  ε = ε_eff
  μ_G_avg = U_s_G_in * D_p * ρ_G_avg / Re_p
  valid = (
    (Re_p >= 2) & (Re_p <= 150) # comiti1989, p. 3
    & (D/D_p >= 10)             # comiti1989, validated for D/D_p = 12 and 53 (Table 2 and D = 60 mm)
    & (phi_s >= 0.95)           # comiti1989, use spherical particles only
  )

  τ = 0.01034 * D_p + 1.38842   # commit1989, tortuosity from Table 4 assuming linear correlation, valid for D_p = 1..5 mm
  α = 6 / D_p                   # commit1989, dynamic specific surface area for spherical particles from p. 5
  M = (1 - 0.0413 * (1 - D_p / D)**2 + 0.0968 * (1 - D_p / D)**2) * τ**2 * ρ_G_avg * α * (1 - ε) / ε**3
  N = 2 * μ_G_avg * τ**2 * α**2 * (1 + 4 / α / D / (1 - ε))**2 * (1 - ε)**2 / ε**3
  ΔP = L * (M * U_s_G_in**2 + N * U_s_G_in)  # Total pressure drop, eq. (17)
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def fand1990(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Fand et al. equation
  https://doi.org/10.1115/1.3242658   1987
  https://doi.org/10.1115/1.2909373   1990
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  M = 1 + 0.666667 * D_p / (1 - ε) / D    # Mehta coeff for wall effects
  Re_w = Re_m / M
  valid = (
    (Re_w >= 3) & (Re_w <= 600)          # fand1990 Fig. 3
    & (D/D_p >= 1.4)                     # fand1990, from Conclusions
    & (phi_s >= 0.95)                    # fand1990
    & (ε >= 0.4) & (ε <= 0.62)           # fand1990 Fig. 2
  )

  data = {
    'A_w': {
      'inf': np.array([172.9, 213.7]),
      'a':   np.array([82.18, 129.7]),
      'p':   np.array([1.125e-4, 3.852e-5]),
      'q':   np.array([-0.003931, -0.003376]),
      'r':   np.array([0.1314, 0.1510])
    },
    'B_w': {
      'inf': np.array([1.871, 1.569]),
      'a':   np.array([1.636, 1.350]),
      'p':   np.array([4.908e-4, 3.688e-4]),
      'q':   np.array([-0.01665, -0.01465]),
      'r':   np.array([0.2925, 0.2646])
    },
    'k': {
      'inf': 5.340, 'a': 0.6545, 'p': 0, 'q': 0, 'r': 0.09034
    }
  }

  def coeffInterpolator(Ylabel):
    idx = (Re_p > 100).astype(int)  # 0 for Re ≤ 100, 1 for Re > 100
    params = data[Ylabel]
    def select(p): return p[idx] if isinstance(p, np.ndarray) else p # If the param is an array, pick based on Re_p. If scalar (like 'k'), use as is.

    inf = select(params['inf'])
    a   = select(params['a'])
    p   = select(params['p'])
    q   = select(params['q'])
    r   = select(params['r'])

    f = p * (D/D_p)**3 + q * (D/D_p)**2 + r * (D/D_p)

    return inf - a * np.exp(-f)

  k = coeffInterpolator('k')              # coeffs for Darcy (laminar) flow
  A = np.where(Re_p < 100, 182, 225)      # coeffs for unbounded packings fand1987
  B = np.where(Re_p < 100, 1.92, 1.61)    # coeffs for unbounded packings fand1987
  A_w = coeffInterpolator('A_w')          # coeffs for bounded packings fand1990
  B_w = coeffInterpolator('B_w')          # coeffs for bounded packings fand1990
  A_eff = np.where(D/D_p < 40, A_w, A)
  B_eff = np.where(D/D_p < 40, B_w, B)

  F_s = np.where(Re_p < 3, 36 * k * M**2 / Re_m, A_eff * M**2 / Re_m + B_eff * M)  # laminar or non-laminar flow
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def foumeny1993(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Foumeny et al. equation
  https://doi.org/10.1016/0017-9310(93)80028-S
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m > 5) & (Re_m < 8500)           # foumeny1993, Fig. 1
    & (D/D_p >= 3)                       # foumeny1993, Fig. 2
    & (phi_s >= 0.95)                    # foumeny1993
    & (ε >= 0.386) & (ε <= 0.467)        # foumeny1993, from Eq. (6)
  )

  B = D / D_p / (0.336 * D / D_p + 2.28)
  F_s = 130 / Re_m + B
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def lee1994(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Lee & Ogawa equation
  https://doi.org/10.1252/jcej.27.691
  """

  D_p = D_p_eff
  ε = ε_eff
  valid = (
    (Re_p > 1) & (Re_p < 300_000)        # lee1994 Figs. 1 and 2
    & (phi_s >= 0.95)                    # lee1994
  )

  n = 0.352 + 0.1 * ε + 0.275 * ε**2    # eq. (9)
  F_s = 6.25 * (29.32 / Re_p + 1.56 / Re_p**n + 0.1) * (1 - ε)     # eq. (11)
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def liu1994(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Liu et al. equation
  https://doi.org/10.1016/0009-2509(94)00168-5
  """

  D_p = D_p_eff
  ε = ε_eff
  μ_G_avg = U_s_G_in * D_p * ρ_G_avg / Re_p
  Re_mod = D_p * ρ_G_avg * U_s_G_in / μ_G_avg * (1 + (1 - ε**0.5)**0.5) / (1 - ε) / ε**0.16667 # eq. (42)
  valid = (
    (Re_mod <= 6000)                     # liu1994, p. 16
    & (D/D_p >= 1.3333)                  # liu1994 from Abstract d_s = D_p/D < 0.75
    & (phi_s >= 0.95)                    # liu1994, p. 15
    & (ε >= 0.36) & (ε <= 0.94)          # liu1994 from Tables 1 and 2 (ΔP correlation was tested on foams to mimic high voidage)
  )

  A = 1 + np.pi * D_p / 6 / (1 - ε) / D
  B = 1 - np.pi**2 * D_p / 24 / D * (1 - 0.5 * D_p / D)
  d_s = D_p / D
  F_v = 85.2 * A**2 + 0.69 * B * Re_mod**3 / (256 + Re_mod**2)
  ΔP = F_v * (μ_G_avg * U_s_G_in * L * (1 - ε)**2) / (d_s**2 * ε**(11/3)) # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def hayes1995(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Hayes et al. equation
  https://doi.org/10.1007/BF01064677
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m >= 0.1) & (Re_m <= 100_000)    # hayes1995 Figs. 2 and 3
    & (D/D_p >= 2)                       # hayes1995 Fig. 5
    & (phi_s >= 0.95)                    # hayes1995
    & (ε >= 0.38) & (ε <= 0.43)          # hayes1995 Table 1
  )
  A = 850
  B = 11.6
  C = 1.3
  T = 1.078147 * ε**3 - 1.46903182 * ε**2 + 0.970189 * ε + 0.248662 # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  F_p = (A + ((B * (3 * T - 1)) * Re_p / (T * (1 - ε) * (1 - T))))**0.5 / T * (1 - ε)**2 / ε**2 + C * Re_p / 2 * T / (3 * T - 1) * (1 - ε) / ε
  ΔP = L * F_p * ρ_G_avg * U_s_G_in**2 / D_p  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def avontuur1996(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Avontuur & Geldart equation
  P.P.C. Avontuur, D. Geldart. A quality assessment of the Ergun equation. *Chem. Eng.* (1996), pp. 994-996.

  !! orig paper wasnt found
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m <= 10_000)                   # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
    #& (D/D_p > 30)                   # assumed
    & (phi_s >= 0.95)                # assumed
    #& (ε >= 0.36) & (ε <= 0.4)       # assumed
  )

  F_s = 141 / Re_m + 1.52  # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def critoph1996 (D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Critoph & Thorpe equation
  https://doi.org/10.1016/1359-4311(95)00023-2
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p >= 30) & (Re_p <= 200)    # critoph1996, Fig. 7 and low-Re statement
    & (L/D >= 2.9)                  # critoph1996, from vessel D=125 mm and bed length ~360 mm
    & (D/D_p >= 50)                 # critoph1996, D=125 mm, d=2.5 mm → D/D_p=50
    & (ε >= 0.36) & (ε <= 0.4)      # critoph1996, measured void fraction ε=0.373 (p. 2)
  )

  F_s = 317 / Re_m + 3.15           # eq. (1)
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def oneill1997(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using O'Neill & Benyahia equation
  K. O'Neill, F. Benyahia. Packed bed systems: an insight into more flexible design. The 1997 IChemE Research Event/The Jubilee Research Event (1997), pp. 1253-1256

  !! orig paper wasnt found
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (D/D_p > 5)                       # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
    & (phi_s >= 0.95)                 # assumed
    #& (ε >= 0.36) & (ε <= 0.4)       # assumed
  )

  A = 521.26 - 22581.24 / (D / D_p)**2
  B = 1.12 + 4.2 / D * D_p
  F_s = A / Re_m + B  # erdim2015revisit https://doi.org/10.1016/j.powtec.2015.06.017
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def raichura1999(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Raichura equation
  https://doi.org/10.1080/089161599269627
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p >= 30) & (Re_p <= 1700)        # raichura1999, Abstract
    & (L/D > 3)                          # raichura1999, deduced from end effects extending ~3 Dp (p.10)
    & (D/D_p >= 5)                       # raichura1999, Abstract
    & (phi_s >= 0.95)                    # raichura1999
    & (ε >= 0.38) & (ε <= 0.43)          # raichura1999 fig. 3
  )

  A = 103 * (ε / (1 - ε))**2 * (6 - 6 * ε + 80 * D_p / D)
  B = 2.8 * ε / (1 - ε) * (1 - 1.82 * D_p / D)**2
  F_s = A / Re_m + B
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def eisfeld2001(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Eisfeld & Schnitzlein equation
  https://doi.org/10.1016/S0009-2509(00)00533-9
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p >= 0.01) & (Re_p <= 17635)     # eisfeld2001, from Conclusions
    & (D/D_p >= 1.624)                   # eisfeld2001, from Conclusions
    & (phi_s >= 0.8)                     # eisfeld2001, correlation for cylinders and spheres
    & (ε > 0.33) & (ε < 0.882)           # eisfeld2001, from Conclusions
  )

  K1 = np.where(phi_s >= 0.95, 154, 155)
  k1 = np.where(phi_s >= 0.95, 1.15, 1.42)
  k2 = np.where(phi_s >= 0.95, 0.87, 0.83)
  A_w = 1 + 0.66667 / (D / D_p) / (1 - ε)
  B_w = (k1 / (D/D_p)**2 + k2)**2
  F_s = K1 * A_w**2 / Re_m + A_w / B_w
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def yu2002(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Yu et al. equation
  https://doi.org/10.1016/S1359-4311(01)00116-8
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p >= 750) & (Re_p <= 2500)   # yu2002
    & (D/D_p >= 30)                  # yu2002 reports D/D_p = 30
    & (phi_s >= 0.95)                # yu2002
    & (ε >= 0.36) & (ε <= 0.38)      # yu2002 reports 0.364 ≤ ε ≤ 0.379, Sec.2.1
  )

  F_s = 203 / Re_m + 1.95
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def felice2004(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, μ_G_avg, **kwargs):
  """
  Calculate pressure drop across packed bed using Di Felice & Gibilaro equation
  https://doi.org/10.1016/j.ces.2004.03.030
  """

  D_p = D_p_eff
  ε = ε_eff
  μ_G_avg = ρ_G_avg * U_s_G_in * D_p / Re_p
  valid = (
    (D/D_p >= 5)                     # felice2004, Fig. 3
    & (phi_s >= 0.95)                # felice2004
    & (ε >= 0.4) & (ε <= 0.5)        # felice2004, from Eq. (14) and D/D_p ≥ 5
  )

  A = 150 * μ_G_avg * (1 - ε)**2 / D_p**2 / ε**3
  B = 1.75 * ρ_G_avg * (1 - ε) / D_p / ε**3
  U_b = U_s_G_in / (2.06 - 1.06 * (D/D_p - 1)**2 / D**2 * D_p**2)
  ΔP = L * (A * U_b + B * U_b**2)
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def nemec2005(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Nemec & Levec equation
  https://doi.org/10.1016/j.ces.2005.05.068
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m >= 1) & (Re_m <= 1000)     # nemec2005, Conclusions
    & (D/D_p >= 10)                  # nemec2005, Sec. 3.2
    & (phi_s >= 0.95)                # nemec2005, A and B coeffs. for spheres only
    & (ε >= 0.35) & (ε <= 0.55)      # nemec2005, Conclusions
  )

  A = 150 / phi_s**1.5
  B = 1.75 / phi_s**1.3333
  F_s = A / Re_m + B
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def montillet2007(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Montillet et al. equation
  https://doi.org/10.1016/j.cep.2006.07.002
  """

  D_p = D_p_eff
  ε = ε_eff
  valid = (
    (Re_p >= 10) & (Re_p <= 2500)        # montillet2007, from Conclusion
    & (D/D_p >= 3.8)                     # montillet2007, from Conclusion
    & (phi_s >= 0.95)                    # montillet2007
    & (ε >= 0.356) & (ε <= 0.452)        # montillet2007, Tables 3 and 4
  )

  a = np.where(ε < 0.4, 0.061, 0.05)
  b = np.where(D/D_p < 50, (D/D_p)**0.2, 2.2)
  F_p = a * b * (1000 / Re_p + 60 / Re_p**0.5 + 12)
  ΔP = L * F_p * ρ_G_avg * U_s_G_in**2 / D_p  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def carpinlioglu2008(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Çarpinlioğlu & Özahi equation
  https://doi.org/10.1016/j.powtec.2008.01.027
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m >= 675) & (Re_m <= 7772)          # carpinliouglu2008simplified, Abstract & Conclusion
    # & (L/D >= 0.24) & (L/D <= 1.46)       # carpinliouglu2008simplified, Section 2, p.2 (verified range)
    & (D/D_p >= 5.72)                       # ccarpinliouglu2008simplified, Table 2
    & (phi_s >= 0.55)                       # ccarpinliouglu2008simplified, Conclusion
    & (ε >= 0.36) & (ε <= 0.56)             # ccarpinliouglu2008simplified, Abstract & Conclusion
  )

  ΔP = 70 * ρ_G_avg * U_s_G_in**2 / (Re_m * D_p / L * ε**7)**0.4733  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def ozahi2008(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Özahi et al. equation
  https://doi.org/10.1163/156855208X314985
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m >= 708) & (Re_m <= 7772)       # ozahi2008, p. 5
    # & (L/D >= 0.24) & (L/D <= 1.46)      # ozahi2008, p. 5 (verified range)
    & (D/D_p >= 5.72)                    # ozahi2008, p. 5
    & (phi_s >= 0.55)                    # ozahi2008, Abstract
    & (ε >= 0.36) & (ε <= 0.56)          # ozahi2008, p. 5
  )

  F_s = 160 / Re_m + 1.61 * phi_s
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def cheng2011(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Cheng equation
  https://doi.org/10.1016/j.powtec.2011.03.026
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p >= 2) & (Re_p <= 5550)        # cheng2011, p. 5
    & (D/D_p >= 1.1)                    # cheng2011, Conclusions
    & (phi_s >= 0.95)                   # cheng2011
    & (ε >= 0.3) & (ε <= 0.7)           # cheng2011 fig. 3
  )

  A = 185 + 17 * ε / (1 - ε) * (D / (D - D_p))**2
  B = 1.3 * ((1 - ε) / ε)**0.3333 + 0.03 * (D / (D - D_p))**2
  F_s = A / Re_m + B
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def harrison2013(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Harrison et al. equation
  https://doi.org/10.1002%2Faic.14034
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p > 0.32) & (Re_p < 7700)        # harrison2013, Summary
    & (D/D_p > 8.3)                      # harrison2013, Summary
    & (phi_s >= 0.95)                    # assumed
    & (ε > 0.33) & (ε < 0.88)            # harrison2013, Summary
  )

  A = 119.8 * (1 + np.pi * D_p / 6 / (1 - ε) / D)**2
  B = 4.63 * (1 - np.pi**2 * D_p / 24 / D * (1 - D_p / 2 / D))
  F_s = A / Re_m + B / Re_m**0.166667
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def erdim2015(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Erdim equation
  https://doi.org/10.1016/j.powtec.2015.06.017
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m > 2) & (Re_m < 3600)         # erdim2015, Conclusions 
    & (D/D_p > 4)                      # erdim2015, Conclusions 
    & (phi_s >= 0.95)                  # erdim2015
    & (ε > 0.37) & (ε < 0.47)          # erdim2015, Conclusions 
  )

  F_s = 160 / Re_m + 2.81 / Re_m**0.096
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def guo2017(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Guo et al. equation
  https://doi.org/10.1016/j.ces.2017.08.022
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (D/D_p >= 1) & (D/D_p <= 2)         # guo2017, Conclusions
    & (phi_s >= 0.95)                   # guo2017
    & (ε > 0.365) & (ε < 0.682)         # guo2017 Table 1
  )

  A = 1004 / (D / D_p)**9.69 + 57.6 * D / D_p
  B = 1964 * D_p / D + 502.7 * D / D_p - 1984
  C = -3.183 * D_p / D - 1.785 * D / D_p + 5.241
  F_s = A / Re_m + B / Re_m**C
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def seckendorff2020(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Seckendorff et al. equation
  https://doi.org/10.1016/j.ces.2020.115644
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_p >= 10) & (Re_p <= 3000)       # seckendorff2020, p. 13
    & (D/D_p >= 4)                      # seckendorff2020, Fig. 15a
    & (phi_s >= 0.86) & (phi_s <= 0.89) # seckendorff2020, For equilateral cylinders ϕ_s = 0.874
    & (ε >= 0.32) & (ε <= 0.45)         # seckendorff2020 Fig. 13a
  )

  F_s = 65.7 / Re_m + 16.25 / Re_m**0.343
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def cheng2021(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Cheng et al. equation
  https://doi.org/10.3390/en14040872
  """

  D_p = D_p_eff
  ε = ε_eff
  valid = (
    (Re_p >= 200) & (Re_p <= 6400)            # cheng2021, Fig. 5
    & (L/D >= 1.86)                           # cheng2021, H=800 mm, D=430 mm → L/D=1.86 (single value, narrow range assumed)
    & (D/D_p >= 13.8)                         # cheng2021, from Table 1 (D=430 mm, d_p values)
    & (phi_s >= 0.69) & (phi_s <= 0.89)       # cheng2021, from Table 1
    & (ε >= 0.52) & (ε <= 0.54)               # cheng2021, from p. 6 (ε=0.53 for d=35 mm; other sizes not reported)
  )

  k = np.where(U_s_G_in < 1.15, 395.2, 17.2)
  a1 = np.where(U_s_G_in < 1.15, -0.47, -0.19)
  a2 = np.where(U_s_G_in < 1.15, -0.5, -0.15)
  F_p = k * Re_p**a1 * (D/D_p)**a2
  ΔP = L * F_p * ρ_G_avg * U_s_G_in**2 / D_p  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def reger2023(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Reger et al. equation
  https://doi.org/10.1016/j.nucengdes.2022.112123
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m >= 100) & (Re_m <= 10000)     # reger2023improved, Conclusion
    & (D/D_p >= 4.4)                    # reger2023improved, Table 1
    & (phi_s >= 0.95)                   # reger2023improved, Table 2
    & (ε >= 0.2) & (ε <= 0.9)           # reger2023improved, Conclusion
  )

  f_ε = 253.9 * ε**4 - 499.3 * ε**3 + 364.7 * ε**2 - 115.6 * ε + 14.21
  F_s = 160 / Re_m + 3 * f_ε / Re_m**0.1
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def dixon2023(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Dixon equation
  https://doi.org/10.1002/aic.18035
  https://www.researchgate.net/publication/367395881_General_correlation_for_pressure_drop_through_randomly-packed_beds_of_spheres_with_negligible_wall_effects
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m >= 0.01) & (Re_m <= 500_000)    # dixon2023, Conclusion
    & (D/D_p >= 5)                        # dixon2023, using eq. (10)
    & (phi_s >= 0.95)                     # dixon2023
  )

  A = 160 * (1 + 2 * 0.5459 / (3 * (1 - ε) * D / D_p ))**2 # Dixon correction that sort of works for low D/Dp (5..15), Eq. (10)
  F_s = A / Re_m + (0.922 + 16 / Re_m**0.46) * Re_m / (Re_m + 52)
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p * (1 - ε) / ε**3  # Total pressure drop
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def wu2025(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Wu & Hibiki equation
  https://doi.org/10.1016/j.ijheatmasstransfer.2024.126620
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  D_H = D_p * ε / 6 / (1 - ε)               # Eq. (6)
  valid = (
    (Re_m >= 0.173) & (Re_m <= 502_000)     # wu2025, Conclusions
    & (L/D_H >= 63.1) & (L/D_H <= 17084)    # wu2025, Conclusions
    & (D_H/D >= 0.0019) & (D_H/D <= 0.0285) # wu2025, Conclusions
    & (phi_s >= 0.95)                       # wu2025, Abstract 
    & (ε >= 0.33) & (ε <= 0.651)            # wu2025, Conclusions 
  )

  F_k_lam = 0.5 + 158 / Re_m**0.8 / (L / D_H)**0.05
  F_k_turb = 0.57 + 2.07 / Re_m**0.12
  F_k = (F_k_lam**3 + F_k_turb**3)**(1/3)
  ΔP = L * F_k * ρ_G_avg * U_s_G_in**2 / D_H  # Total pressure drop, eq. (4)
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }

############################################################################################
def xu2026(D, L, D_p_eff, ε_eff, Re_p, phi_s, ρ_G_avg, U_s_G_in, **kwargs):
  """
  Calculate pressure drop across packed bed using Xu et al. equation
  https://doi.org/10.1016/j.powtec.2025.121784
  """

  D_p = D_p_eff
  ε = ε_eff
  Re_m = Re_p / (1 - ε)
  valid = (
    (Re_m >= 350) & (Re_m <= 4980)         # xu2026, Sec. 4.2, Eq. (13) applicability statement
    & (L/D >= 1.8)                         # xu2026, deduced from L=600-1000 mm, D=219-325 mm (Sec. 2.1)
    & (D/D_p >= 13)                        # xu2026, deduced from D=219-325 mm, d_sp=9.15-14.51 mm (Table 2, Sec. 2.1)
    & (phi_s >= 0.578) & (phi_s <= 0.846)  # xu2026, abstract and conclusion
    & (ε >= 0.57) & (ε <= 0.64)            # xu2026, abstract
  )

  k1 = 902.27 - 43.33 * np.exp(D_p / 22.28)
  k2 = 8.7 - 0.3 * np.exp(D_p / 14.8)
  F_s = k1 / Re_m / phi_s + k2
  ΔP = L * F_s * ρ_G_avg * U_s_G_in**2 / D_p / phi_s * (1 - ε) / ε**3  # Total pressure drop, Eq. (13)
  ΔP_array = np.where(valid, ΔP, np.nan)

  return { 'ΔP': ΔP_array }