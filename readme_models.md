# Model Types

**analytical** — pure math without any awkward coeffs\
**semi-empirical** — elegant generalized correlations but with awkward coeffs\
**empirical** — meaningless correlations with tight validity ranges

***

# Pressure Drop Models

Vector functions that accept either single floats or NumPy arrays.

- ϕ\_s ≥ 0.95 → spherical
- Models not designed for irrigated packings use `D_p_eff` and `ε_eff` values to mimic liquid holdup
- Upper limits for $D/D\_p$ and $L/D$ are omitted except special cases when correlations are for narrow/shallow beds

Both laminar (`# Laminar Stress Test`) and turbulent (`# Turbulent Stress Test`) regimes are handled via A-type rules using local Spearman correlation with expected positive correlation (∂ΔP/∂X > 0).

**The following models were adopted without confirmation by original texts:**

1. morcom1946 (src: [erdim2015revisit](https://doi.org/10.1016/j.powtec.2015.06.017))
2. leva1949 (src: [erdim2015revisit](https://doi.org/10.1016/j.powtec.2015.06.017))
3. fahien1961 (src: [erdim2015revisit](https://doi.org/10.1016/j.powtec.2015.06.017))
4. handley1968 (src: [erdim2015revisit](https://doi.org/10.1016/j.powtec.2015.06.017))
5. kuo1978(src: [erdim2015revisit](https://doi.org/10.1016/j.powtec.2015.06.017) and [DTIC](https://apps.dtic.mil/sti/tr/pdf/ADA091300.pdf))
6. paterson1986 (src: [erdim2015revisit](https://doi.org/10.1016/j.powtec.2015.06.017))
7. watanabe1989 (src: [erdim2015revisit](https://doi.org/10.1016/j.powtec.2015.06.017))
8. avontuur1996 (src: [erdim2015revisit](https://doi.org/10.1016/j.powtec.2015.06.017))
9. oneill1997 (src: [erdim2015revisit](https://doi.org/10.1016/j.powtec.2015.06.017))

**The following models were adopted for packed beds:**

1. kurten1966 (adopted by [ozahi2008](https://doi.org/10.1163/156855208X314985))

***

# Void Fraction Models

- `benyahia2005` — Benyahia & O'Neill equation for mean void fraction of a packed bed
- `ribeiro2010` — Ribeiro et al. equation for mean void fraction of a packed bed

Zou & Yu (1995) void fraction correlation for dense and loose packings of spheres is also referenced.

***

# Mass Transfer Models

Vector functions for mass transfer calculations in packed bed condensers.

- `massHTU` — calculate mass Height of Transfer Unit
- `outputGasComposition` — calculate molar fractions of output gas mixture
- `calculate_G_cond` — calculate condensate mass flow rate from gas humidity change

**Assumptions:** gas exits saturated (RH=100%) at water inlet temperature; dry gas composition ratio (CO₂/N₂/O₂) assumed constant.

***

# Heat Transfer

Not yet implemented (placeholder).

***

# Data Mapping & Derived Parameters

**Input params** — packed bed geometry (D, L, D\_p, ϕ\_s, ε), flue gas composition (M\_G\_in\_CO₂, M\_G\_in\_H₂O, M\_G\_in\_O₂), gas velocity & temperature (U\_s\_G\_in, T\_G\_in), liquid parameters (Ratio\_L\_\_G, T\_wtr\_in).

**Output params** — ΔP, T\_G\_out, T\_L\_out, M\_G\_out\_CO₂, M\_G\_out\_H₂O, M\_G\_out\_O₂.

**Derived (Data Map):** T\_avg, ρ\_G\_avg, μ\_G\_avg, ρ\_L\_avg, μ\_L\_avg, G\_cond, h\_L (liquid holdup), ε\_eff (wet void fraction), D\_p\_eff (wet particle diameter), Re\_p, Re\_m.
