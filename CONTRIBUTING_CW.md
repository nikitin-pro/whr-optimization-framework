# Contributing to the Cross-Validation Framework

Thank you for your interest in contributing to the Cross-Validation Framework. This document outlines the guidelines for contributing new models to the framework.

## Model Contribution Requirements

Each model contribution must include the following information. Please follow this template when adding a new model to `src/models_pressure_drop.py` and registering it in `src/config.py`.

### 1. Model Type

Specify the type of the model based on its methodology:

- **analytical**: Pure math without any empirical coefficients
- **semi-empirical**: Elegant generalized correlations with empirical coefficients
- **empirical**: Correlations with tight validity ranges and limited physical basis

### 2. List of Input Parameters

Provide a list of input parameters in the optimization context. Common input parameters include:

- `Re_p` — Particle Reynolds number
- `L` —Packed bed length
- `D` — Packed bed equivalent diameter
- `phi_s` — Particle sphericity
- `ε` — Void fraction
- `h_L` - Liquid holdup (for two-phase models)

> Note that $Re_p=\frac{\rho U_sD_p}{\mu}$. Don't specify $Re_p$ along with its components.

### 3. List of Output Parameters

Provide a list of output parameters that the model calculates. Common output parameters include:

- `ΔP` - Pressure drop

### 4. Model Function

Implement the model function following the template below (based on xu2026):

```python
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
```

### 5. Source Publication DOI

Provide the DOI (Digital Object Identifier) of the original publication that introduced the model. This is essential for proper attribution and for users to access the original source material.

## Registration in config.py

After implementing your model function, register it in `src/config.py` by adding an entry to the `src_models` list:

```python
{'name':'your_model_name', 'src_type':'semi-empirical', 'inContext':['Re_p', 'L'], 'outContext':['ΔP']},
```

Also add the function name to the `__all__` list in `src/models_pressure_drop.py`.

## Quality Checklist

Before submitting a contribution, please ensure:

- [ ] Model function follows the established template and conventions.
- [ ] All validity ranges are clearly documented with references to the original publication.
- [ ] DOI is provided and verified.
- [ ] Model is registered in `src/config.py` with correct `inContext` and `outContext` parameters.
- [ ] Function name is added to `__all__` in `src/models_pressure_drop.py`.
- [ ] Any special parameters (e.g., `h_L` for liquid holdup) are properly handled.

## Testing

Ensure your model works correctly by:

1. Verifying the function runs without errors.
2. Checking that validity conditions properly mask invalid inputs.
3. Comparing outputs against known values from the original publication where possible.
