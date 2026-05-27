# Optimization Framework

## Concept

### 1. Dataset Generation

1.1. Generate input design space (X) using Sobol with n dimensions where n is the number of inputs.\
1.2. Create a dynamic data map that will contain all inputs X\_i and also derived inputs (e.g., Us, Dp, ρ, μ → Re).\
1.3. Generate extended output dataframe Y\_ext (len(Y\_ext) = len(X), column\_names = model\_names) by calculating each model over the whole design space X.\
1.4. Generate validity dataframe V\_ext (shaped as Y\_ext) and populate it with binary validities (1 if inside validity range, NaN otherwise), claimed by the models.

### 2. Apply A-Rules (1 input X\_i → 1 output Y\_i)

E.g., L → ΔP or Re=f(Us,Dp,ρ,μ) → ΔP. For each model:\
2.1. Determine design subspace from validity ranges claimed by model.\
2.2. Select anchor points (about 10% of subspace; distributed).\
2.3. Calculate validity score of the model in these anchor points for each a-rule:\
2.3.1. For each DP (row of the design subspace), variate the input (e.g., L or Re) from a-rule (10 samples uniformly distributed over the subspace in the given dimension, e.g. 0.1\<L<10 m or 10\<Re<6000) → calculate model output (e.g., ΔP) trend (gradient), e.g., ∂ΔP/∂L or ∂ΔP/∂Re.\
2.3.2. Compare this gradient with expected correlation sign (positive/negative) from the a-rule → validity score of the a-rule.\
2.4. Consolidate individual validity scores from a-rules for each model.\
2.5. Map consolidated validity scores onto the whole design space V\_ext using classifier.

### 3. Apply B-Rules (Conservation Laws)

Multiple inputs and outputs → residual, e.g., T\_in, T\_out, Gm, Pr → Q\_in - Q\_out = residual.\
3.1. For each model, calculate validity score of the model (convert from residual) in every DP (row of the whole design space) for each b-rule.\
3.2. Update V\_ext by consolidating with validity scores from b-rules.

### 4. Initial Consensus

4.1. Calculate consensus outputs (Y) using medians of individual results from Y\_ext weighted with validity scores from V\_ext.\
4.2. Consolidate validity scores from V\_ext into dataframe V shaped as Y.

### 5. Check Alignment of the Models

For each DP, compare individual model results from Y\_ext with consensus from Y → put alignment scores into A\_ext shaped as Y\_ext.

### 6. Final Consensus

6.1. Update Y using Y\_ext, updated V\_ext, and A\_ext. Discard models with low alignment score if there are >10 overlapping models for a given DP. Empirical models should be treated superior to semi-empirical.\
6.2. Update V using V\_ext and A\_ext (multiply V\_ext and A\_ext; empirical models should be treated superior to semi-empirical). Exclude discarded models

### 7. Train Surrogate Model

Train surrogate model on X, Y, and V.
