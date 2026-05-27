import numpy as np
import pandas as pd
from scipy.stats import qmc

from config import input_params, output_params, src_models
from models_pressure_drop import *
from models_mass_transfer import *
from models_heat_transfer import *
from models_misc import *


def run():
    print("\033[1m" + "=" * 60)
    print("SECTION 1: Dataset Generation")
    print("=" * 60 + "\033[0m")

    # 1.1 Input Space (X) — Sobol Sequence
    print("\nGenerating input space using Sobol sequence...")
    sampler = qmc.Sobol(d=len(input_params), scramble=True)
    sample = sampler.random_base2(m=18)
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
    print('Input parameters:', list(X.columns))

    # 1.2 Data Map
    print("\nBuilding data map...")
    data_map = {col: X[col].values for col in X.columns}
    data_map['T_avg'] = ((X['T_G_in'] + X['T_wtr_in']) / 2).values
    data_map['ρ_G_avg'], data_map['μ_G_avg'] = gas_props(
        M_CO2=X['M_G_in_CO2'], M_H2O=X['M_G_in_H2O'],
        M_O2=X['M_G_in_O2'], T=data_map['T_avg'], props=['ρ', 'μ']
    )
    data_map['ρ_L_avg'], data_map['μ_L_avg'] = water_props(
        T=data_map['T_avg'], props=['ρ', 'μ']
    )
    data_map['G_cond'] = calculate_G_cond(**{   # estimation of condensate mass flow rate (excluding rinse water)
        col: X[col].values for col in X.columns
        if col in ['M_G_in_H2O', 'M_G_in_CO2', 'M_G_in_O2', 'U_s_G_in', 'T_G_in', 'T_wtr_in', 'D']
    })
    data_map['h_L'] = liquid_holup(**data_map)
    data_map['ε_eff'] = X['ε'] - data_map['h_L']  # effective (wet bed) void fraction (reduced by volume occupied by the liquid)
    data_map['D_p_eff'] = X['D_p'] * ((1 - data_map['ε_eff']) / (1 - X['ε'])) ** 0.3333  # effective (wet bed) particle diameter (enlarged by liquid volume) https://doi.org/10.1016/0950-4214(89)80016-7
    data_map['Re_p'] = data_map['ρ_G_avg'] * X['D_p'] * X['U_s_G_in'] / data_map['μ_G_avg']
    data_map['Re_m'] = data_map['Re_p'] / (1 - X['ε'])

    print("\033[32mData map initialized\033[0m\n")
    print('Data map keys:', list(data_map.keys()))
    
    print("\nRe_p range: [{:.2e}, {:.2e}]".format(data_map['Re_p'].min(), data_map['Re_p'].max()))
    print("Re_m range: [{:.2e}, {:.2e}]".format(data_map['Re_m'].min(), data_map['Re_m'].max()))
    print("L/D range: [{:.2f}, {:.2f}]".format((X['L'] / X['D']).min(), (X['L'] / X['D']).max()))
    print("D/D_p range: [{:.2f}, {:.2f}]".format((X['D'] / X['D_p']).min(), (X['D'] / X['D_p']).max()))

    # Flooding Check
    is_valid = (data_map['h_L'] / X['ε'].values) <= 0.5
    data_valid = {k: v[is_valid] for k, v in data_map.items()}
    print(f"Valid design points: {np.sum(is_valid)} / {len(X)}")

    # 1.3 Output Space (Y_ext)
    print("\nCalculating source model predictions...")
    import models_pressure_drop, models_mass_transfer, models_heat_transfer, models_misc
    _model_funcs = {}
    for _mod in [models_pressure_drop, models_mass_transfer, models_heat_transfer, models_misc]:
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
    for col in Y_ext.columns.tolist():
        model_samples = Y_ext[col].count()
        print(f"  {col:<25} {model_samples:<6} ({(model_samples / total_samples * 100):.2f}%)")

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