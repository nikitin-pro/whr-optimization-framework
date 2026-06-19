import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import qmc

from .config import input_params, output_params, src_models
from .models_pressure_drop import *
from .models_misc import *


def run():
    print("\033[1m" + "=" * 60)
    print("SECTION 1: Dataset Generation")
    print("=" * 60 + "\033[0m")

    # 1.1 Input Space (X) — Sobol Sequence
    print("\nGenerating input space using Sobol sequence...")
    sampler = qmc.Sobol(d=len(input_params), scramble=True)
    sample = sampler.random_base2(m=18)
    # Convert open intervals to closed intervals [0, 1] for every dimension
    # axis=0 applies the min-max scaling column-by-column
    sample = (sample - sample.min(axis=0)) / (sample.max(axis=0) - sample.min(axis=0))
    X = pd.DataFrame()
    for i, p in enumerate(input_params):
        xi = sample[:, i]

        if p['type'] == 'independent':
            low, high = p['range']
            if p['logScaling']:
                log_min, log_max = np.log10(low), np.log10(high)
                X[p['name']] = 10**(log_min + xi * (log_max - log_min))
            else:
                X[p['name']] = low + xi * (high - low)

        elif p['type'] == 'dependent':
            base_val = X[p['base']]
            low, high = p['ratio_range']
            if p['logScaling']:
                log_min = np.log10(base_val * low)
                log_max = np.log10(base_val * high)
                X[p['name']] = 10**(log_min + xi * (log_max - log_min))
            else:
                X[p['name']] = (base_val * low) + xi * (base_val * high - base_val * low)

    print("\033[32mInput Space initialized\033[0m\n")
    print('Samples:', len(X))
    print('Input parameters:', [c.replace('\u03b5', 'e') if isinstance(c, str) else c for c in X.columns])

    # 1.2 Data Map
    print("\nBuilding data map...")
    data_map = {col: X[col].values for col in X.columns}
    data_map['ρ_G_avg'], data_map['μ_G_avg'] = gas_props(
        M_CO2=np.full(len(X), 0.0004), M_H2O=np.full(len(X), 0.0115),
        M_O2=np.full(len(X), 0.209), T=np.full(len(X), 20), props=['ρ', 'μ']
    ) # simplified air at 20 ℃ and 1 atm
    data_map['h_L'] = np.zeros(len(X)) # placeholder for liquid holdup
    data_map['ε_eff'] = X['ε'] # placeholder for effective (wet bed) void fraction
    data_map['D_p_eff'] = X['D_p']  # placeholder for effective (wet bed) particle diameter
    data_map['Re_p'] = data_map['ρ_G_avg'] * X['D_p'] * X['U_s_G_in'] / data_map['μ_G_avg']
    data_map['Re_m'] = data_map['Re_p'] / (1 - X['ε'])

    print("\033[32mData map initialized\033[0m\n")
    print('Data map keys:', list(data_map.keys()))
    
    print("\nRe_p range: [{:.2e}, {:.2e}]".format(data_map['Re_p'].min(), data_map['Re_p'].max()))
    print("Re_m range: [{:.2e}, {:.2e}]".format(data_map['Re_m'].min(), data_map['Re_m'].max()))
    print("L/D range: [{:.2f}, {:.2f}]".format((X['L'] / X['D']).min(), (X['L'] / X['D']).max()))
    print("D/D_p range: [{:.2f}, {:.2f}]".format((X['D'] / X['D_p']).min(), (X['D'] / X['D_p']).max()))


    is_valid = np.ones(len(X), dtype=bool) # skip flooding check — all points are valid
    data_valid = {k: v[is_valid] for k, v in data_map.items()}
    print(f"Valid design points: {np.sum(is_valid)} / {len(X)}")

    # 1.3 Output Space (Y_ext)
    print("\nCalculating source model predictions...")
    from . import models_pressure_drop, models_misc
    _model_funcs = {}
    for _mod in [models_pressure_drop, models_misc]:
        for _name in _mod.__all__:
            _model_funcs[_name] = getattr(_mod, _name)

    results_dict = {}
    for model in src_models:
        model_name = model['name']
        model_func = _model_funcs.get(model_name)
        if model_func is None:
            continue
        results = model_func(**data_valid)
        for out_name, out_array in results.items():
            full_col = np.full(len(X), np.nan)
            full_col[is_valid] = out_array
            results_dict[f"{model['name']}__{out_name}"] = full_col

    Y_ext = pd.DataFrame(results_dict, index=X.index)
    print("\033[32mY_ext Dataframe initialized\033[0m\n")
    print("\033[1mOutput parameters:\033[0m")
    total_samples = len(X)

    print(f"\n\033[1mDesign Space Occupancy:\033[0m")
    for out in output_params:
        cols = [col for col in Y_ext.columns if col.endswith(f"__{out['name']}")]
        count = Y_ext[cols].notna().any(axis=1).sum()
        print(f"  Design space coverage for {out['name']:<12} is {(count/total_samples*100):.2f}%")

    # Initialize Validity Dataframes
    output_param_names = [p['name'] for p in output_params]
    col_map = {c: c.split('__')[-1] for c in Y_ext.columns if c.split('__')[-1] in output_param_names}
    V_ext = Y_ext.notna().astype(float).replace(0, np.nan)
    V = V_ext.T.groupby(col_map).max().T

    # 1.4 Design Space Coverage Histogram (Cumulative)
    print("\n\033[1mDesign Space Coverage Histogram:\033[0m")
    for out in output_params:
        out_name = out['name']
        cols = [col for col in V_ext.columns if col.endswith(f"__{out_name}")]
        if not cols:
            continue

        # Count how many models produce a valid prediction for each sample
        overlap_counts = V_ext[cols].notna().sum(axis=1)
        n_total = len(overlap_counts)

        # Cumulative bins: 0, 1+, 2+, ..., 10+
        max_k = min(10, len(cols))
        bins = list(range(max_k + 1))  # 0, 1, 2, ..., 10
        labels = ['0'] + [f'{k}+' for k in range(1, max_k + 1)]

        cumulative_pcts = []
        for k in bins:
            if k == 0:
                pct = (overlap_counts == 0).sum() / n_total * 100
            else:
                pct = (overlap_counts >= k).sum() / n_total * 100
            cumulative_pcts.append(pct)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(bins, cumulative_pcts, width=0.6, edgecolor='black', color='steelblue')
        ax.set_xlabel('Minimum Number of Overlapping Models')
        ax.set_ylabel('Percentage of Design Space [%]')
        ax.set_title(f'Cumulative Design Space Coverage — {out_name}')
        ax.set_xticks(bins)
        ax.set_xticklabels(labels)
        ax.set_xlim(-0.5, max_k + 0.5)

        for x, y in zip(bins, cumulative_pcts):
            ax.text(x, y + 0.3, f'{y:.1f}%', ha='center', va='bottom', fontsize=7)

        plt.tight_layout()
        plt.show()

    return {
        'X': X,
        'data_map': data_map,
        'is_valid': is_valid,
        'Y_ext': Y_ext,
        'V_ext': V_ext,
        'V': V,
    }


if __name__ == '__main__':
    outputs = run()
    print("\033[32mDataset generation complete.\033[0m")
    print(f"X shape: {outputs['X'].shape}")
    print(f"Y_ext shape: {outputs['Y_ext'].shape}")
    print(f"Valid points: {np.sum(outputs['is_valid'])}")