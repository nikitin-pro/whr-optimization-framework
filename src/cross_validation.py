import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import rankdata
import bottleneck as bn
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import re

from .config import input_params, output_params, scaling_lookup, src_models, a_rules, b_rules
from .models_misc import gas_props, water_props, liquid_holup, calculate_G_cond

penalty_threshold = 8 # theta; used in weighted_median function, check_model_coverage function, and model coverage output (run function)

def pretty_label(raw):
    """Format raw model name like 'erdim2015' into 'Erdim (2015)'."""
    import re
    m = re.match(r'^(.*?)(\d{4})$', raw)
    if m:
        return f"{m.group(1).capitalize()} ({m.group(2)})"
    return raw.capitalize()


def color_penalties(val):
    if val == "X":
        return 'color: gray; font-weight: bold'
    try:
        f_val = float(val)
        if f_val > 50:
            return 'color: red; font-weight: bold'
        elif 20 <= f_val <= 50:
            return 'color: #D4AC0D; font-weight: bold'
        else:
            return 'color: black'
    except:
        return 'color: black'


def weighted_median(values_matrix, weights_matrix):
    N, M = values_matrix.shape
    consensus = np.full(N, np.nan)
    confidence = np.zeros(N)
    w_eps = weights_matrix + 1e-6
    for i in range(N):
        v = values_matrix[i]
        w = w_eps[i]
        mask = ~np.isnan(v)
        v_m, w_m = v[mask], w[mask]
        n = len(v_m)
        if n == 0:
            continue
        if n <= 2:
            consensus[i] = np.sum(v_m * w_m) / np.sum(w_m)
        else:
            sorter = np.argsort(v_m)
            v_s, w_s = v_m[sorter], w_m[sorter]
            consensus[i] = v_s[np.searchsorted(np.cumsum(w_s), 0.5 * np.sum(w_s))]
        log_n = np.log(n + 1)
        log_threshold = np.log(penalty_threshold + 1)
        density_penalty = min(log_n / log_threshold, 1.0)
        avg_vld = np.mean(w_m)
        agreement = 1
        if n > 1:
            std_v = np.std(v_m)
            mean_v = np.abs(np.mean(v_m)) + 1e-6
            agreement = np.exp(-(std_v / mean_v))
        confidence[i] = avg_vld * agreement * density_penalty
    return consensus, confidence


def energyBalanceResidual(df_x, df_y):
    q_hot = df_x['L'] * df_y['ΔP']
    q_cold = df_x['L'] * df_x['ε'] * df_y['ΔP']
    return q_hot - q_cold


#deprecated function
def check_rule_overlaps(rules):
    registry = {}
    def is_overlap(r1, r2):
        s1 = r1[0] if r1[0] is not None else float('-inf')
        e1 = r1[1] if r1[1] is not None else float('inf')
        s2 = r2[0] if r2[0] is not None else float('-inf')
        e2 = r2[1] if r2[1] is not None else float('inf')
        return s1 < e2 and s2 < e1
    for rule in rules:
        io_pair = (rule['input'], rule['output'])
        if io_pair not in registry:
            registry[io_pair] = []
        corrs = rule['correlations']
        for i, c1 in enumerate(corrs):
            for c2 in corrs[i + 1:]:
                if is_overlap(c1['range'], c2['range']) and c1['corr'] != c2['corr']:
                    raise ValueError(
                        f"\033[91mConflict WITHIN rule '{rule['name']}': \033[0m"
                        f"Range {c1['range']} ({c1['corr']}) overlaps with {c2['range']} ({c2['corr']})"
                    )
        for entry in corrs:
            new_r, new_c = entry['range'], entry['corr']
            for old_r, old_c, old_name in registry[io_pair]:
                if is_overlap(new_r, old_r) and new_c != old_c:
                    raise ValueError(
                        f"\033[91mConflict BETWEEN '{old_name}' and '{rule['name']}': \033[0m"
                        f"Range {old_r} ({old_c}) overlaps with {new_r} ({new_c}) for {io_pair}"
                    )
            registry[io_pair].append((new_r, new_c, rule['name']))
    print(f"\033[32mSuccess: No conflicting overlaps found.\033[0m\n")


def calculate_a_rule_validity(rule, active_mask, model, data_map, Y_ext, X):
    N_active = np.sum(active_mask)
    scores_active = np.full(N_active, 1.0)
    ctx = {
        'dm': {k: v[active_mask] for k, v in data_map.items()},
        'y': Y_ext[f"{model['name']}__{rule['y']}"].values[active_mask]
    }
    X_test = rule['input'](**ctx)
    Y_test = rule['output'](**ctx)
    ranging_data = data_map[rule['control_param']][active_mask]
    low, high = rule['control_range']
    low, high = (low if low is not None else -np.inf), (high if high is not None else np.inf)
    range_mask = (ranging_data >= low) & (ranging_data <= high)
    if not np.any(range_mask):
        final_scores = np.full(len(X), np.nan)
        final_scores[active_mask] = scores_active
        return final_scores
    context_data = []
    for feat in model.get('inContext', []):
        if scaling_lookup.get(feat, False):
            val = np.log10(np.asarray(data_map[feat])[active_mask] + 1e-6)
        else:
            val = np.asarray(data_map[feat])[active_mask]
        context_data.append(val)
    X_context = np.stack(context_data, axis=1)
    scaler = MinMaxScaler()
    X_scaled_context = scaler.fit_transform(X_context)
    tree = BallTree(X_scaled_context, leaf_size=100)
    k = max(40, int(np.sqrt(N_active)))
    distances, indices = tree.query(X_scaled_context, k=k)

    def batch_spearman_simple(a, b):
        a_ranked = bn.nanrankdata(a, axis=1)
        b_ranked = bn.nanrankdata(b, axis=1)
        mask = np.isfinite(a) & np.isfinite(b)
        a_r = np.where(mask, a_ranked, np.nan)
        b_r = np.where(mask, b_ranked, np.nan)
        mu_a = bn.nanmean(a_r, axis=1)
        mu_b = bn.nanmean(b_r, axis=1)
        a_m = a_r - mu_a[:, np.newaxis]
        b_m = b_r - mu_b[:, np.newaxis]
        num = np.nansum(a_m * b_m, axis=1)
        den = np.sqrt(bn.nansum(a_m**2, axis=1) * bn.nansum(b_m**2, axis=1))
        return np.divide(num, den + 1e-6, out=np.zeros_like(num))

    X_nb = X_test[indices]
    Y_nb = Y_test[indices]
    X_nb_ranked = rankdata(X_nb, axis=1, method='average')
    Y_nb_ranked = rankdata(Y_nb, axis=1, method='average')
    corr_vals_active = batch_spearman_simple(X_nb_ranked, Y_nb_ranked)
    avg_dists = np.mean(distances, axis=1)
    reliability_limit = np.nanmedian(avg_dists) * 2
    unreliable_mask = avg_dists > reliability_limit
    unique_y_per_nb = np.array([len(np.unique(Y_nb[i])) for i in range(len(Y_nb))])
    c_v = corr_vals_active[range_mask]
    local_sensitivity = np.nanmean(np.abs(c_v))
    if local_sensitivity < 0.1 or np.isnan(local_sensitivity):
        scores_active[range_mask] = np.nan
    else:
        u_v = unique_y_per_nb[range_mask]
        local_avg_uniqueness = np.nanmean(u_v) / k
        if local_avg_uniqueness < 0.3:
            scores_active[range_mask] = np.nan
        else:
            if rule['expected_corr'] == 'pos':
                scores_active[range_mask] = np.where(c_v > 0.5, 1.0, np.maximum(0, c_v * 2))
            elif rule['expected_corr'] == 'neg':
                scores_active[range_mask] = np.where(c_v < -0.5, 1.0, np.maximum(0, -c_v * 2))
    scores_active[unreliable_mask] = np.nan
    scores_active[np.isnan(corr_vals_active)] = np.nan
    final_scores = np.full(len(X), np.nan)
    final_scores[active_mask] = scores_active
    return final_scores


def calculate_b_rule_validity(df_X, df_Y, V, Y_ext, V_ext):
    baseline_res = abs(energyBalanceResidual(df_X, df_Y))
    for col in Y_ext.columns:
        p_name = col.split('__')[-1]
        if p_name not in [p['name'] for p in output_params]:
            continue
        p_cols = [c for c in Y_ext.columns if c.endswith(f"__{p_name}")]
        n_active = (V_ext[p_cols] > 0).sum(axis=1)
        test_Y = df_Y.copy()
        test_Y[p_name] = Y_ext[col]
        model_res = abs(energyBalanceResidual(df_X, test_Y))
        res_ratio = np.where(model_res <= baseline_res, 1.0,
                             baseline_res / (model_res + 1e-6))
        discard_mask = (n_active >= 10) & (res_ratio < 0.5)
        V.loc[discard_mask, col] = 0.0
        penalize_mask = (n_active < 10)
        V.loc[penalize_mask, col] *= (res_ratio[penalize_mask] * 0.1)
        V.loc[~discard_mask & ~penalize_mask, col] *= res_ratio[~discard_mask & ~penalize_mask]
    return V


def check_model_alignment(Y_ext, Y, V_ext, threshold_discard=0.5):
    V_aligned = V_ext.copy()
    for col in Y_ext.columns:
        p_name = col.split('__')[-1]
        if p_name not in Y.columns:
            continue
        baseline = Y[p_name].values
        model_val = Y_ext[col].values
        deviation = np.abs(model_val - baseline) / (np.abs(baseline) + 1e-6)
        p_cols = [c for c in Y_ext.columns if c.endswith(f'__{p_name}')]
        n_active = Y_ext[p_cols].notna().sum(axis=1).values
        discard_mask = (n_active >= penalty_threshold) & (deviation > threshold_discard)
        V_aligned.loc[discard_mask, col] = 0.0
        penalize_mask = (n_active < penalty_threshold) & (deviation > threshold_discard)
        V_aligned.loc[penalize_mask, col] *= 0.2
    return V_aligned


def plot_ensemble_validity(V_ext, Y_ext, p_desc, output_dir='.'):
    models_cols = V_ext.columns
    n_models = len(models_cols)
    n_bins = 20
    bin_edges = np.linspace(0, 1, n_bins + 1)
    total_ds_points = len(Y_ext)
    if n_models == 0:
        print(f"No models to plot for {p_desc}")
        return
    sort_idx = V_ext.mean().sort_values(ascending=True).index
    fig, ax = plt.subplots(figsize=(14, 0.4 * n_models + 2))
    for i, col in enumerate(sort_idx):
        model_data_mask = Y_ext[col].notna()
        n_model_data = model_data_mask.sum()
        if n_model_data == 0:
            continue
        valid_scores = V_ext[col].dropna().values
        n_verified = len(valid_scores)
        if n_verified == 0:
            continue
        pct_fail = (np.sum(valid_scores < 0.2) / n_verified * 100)
        pct_high = (np.sum(valid_scores > 0.8) / n_verified * 100)
        pct_mid = (np.sum((valid_scores >= 0.2) & (valid_scores <= 0.8)) / n_verified * 100)
        pct_unval = ((n_model_data - n_verified) / n_model_data) * 100
        counts, _ = np.histogram(valid_scores, bins=bin_edges)
        bar_pcts = (counts / n_model_data) * 100
        left = 0
        for j, pct in enumerate(bar_pcts):
            if pct > 0:
                if bin_edges[j] < 0.2:
                    color = 'salmon'
                elif bin_edges[j] < 0.8:
                    color = 'skyblue'
                else:
                    color = 'lightgreen'
                ax.barh(i, pct, left=left, height=0.7, color=color, edgecolor='black', linewidth=0.2)
            left += pct
        ax.text(101, i, f"{pct_fail:.0f}%", color='salmon', va='center', fontweight='bold', fontsize=10)
        ax.text(106, i, f"{pct_mid:.0f}%", color='skyblue', va='center', fontweight='bold', fontsize=10)
        ax.text(111, i, f"{pct_high:.0f}%", color='lightgreen', va='center', fontweight='bold', fontsize=10)
        ax.text(116, i, f"{pct_unval:.0f}%", color='black', va='center', fontsize=10)
    y_labels = []
    for m in sort_idx:
        ds_pct = (Y_ext[m].notna().sum() / total_ds_points) * 100
        y_labels.append(f"{pretty_label(m.split('__')[0])} ({ds_pct:.2f}% of DS)")
    ax.set_yticks(range(n_models))
    ax.set_yticklabels(y_labels, fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_xlabel('Verified Model Data (% of Claimed Range)', fontsize=12)
    ax.set_title('Physics Validity Distribution per Model', fontsize=14, fontweight='bold', pad=20)
    ax.text(0.5, 1.02, p_desc, fontsize=12, fontstyle='italic', ha='center', va='bottom', transform=ax.transAxes)
    legend_elements = [
        Line2D([0], [0], color='salmon', lw=6, label='Low Validity (<0.2)'),
        Line2D([0], [0], color='skyblue', lw=6, label='Medium Validity (0.2-0.8)'),
        Line2D([0], [0], color='lightgreen', lw=6, label='High Validity (>0.8)'),
        Line2D([0], [0], color='black', lw=6, label='Not Validated (Skipped)')
    ]
    ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)
    # Collect CSV data
    csv_rows = []
    for m in sort_idx:
        model_data_mask = Y_ext[m].notna()
        n_model_data = model_data_mask.sum()
        if n_model_data == 0:
            continue
        valid_scores = V_ext[m].dropna().values
        n_verified = len(valid_scores)
        if n_verified == 0:
            continue
        pct_low = round(np.sum(valid_scores < 0.2) / n_verified * 100, 2)
        pct_med = round(np.sum((valid_scores >= 0.2) & (valid_scores <= 0.8)) / n_verified * 100, 2)
        pct_high = round(np.sum(valid_scores > 0.8) / n_verified * 100, 2)
        label = pretty_label(m.split('__')[0])
        csv_rows.append({'Label': label, 'High': pct_high, 'Med': pct_med, 'Low': pct_low})

    # Save CSV
    if csv_rows:
        csv_filename = f"{output_dir}/validity_{p_desc.replace(' ', '_').replace('℃', 'C')}.csv"
        df_csv = pd.DataFrame(csv_rows)
        df_csv = df_csv.sort_values('High', ascending=False)
        df_csv.to_csv(csv_filename, index=False)
        print(f"CSV saved to {csv_filename}")

    plt.grid(axis='x', linestyle='--', alpha=0.3)
    plt.subplots_adjust(right=0.85)
    plt.savefig(f"{output_dir}/validity_{p_desc.replace(' ', '_').replace('℃', 'C')}.png", dpi=150, bbox_inches='tight')
    plt.close()


def run(ds_outputs):
    print("\033[1m" + "=" * 60)
    print("SECTION 2: Cross-Validation")
    print("=" * 60 + "\033[0m")

    X = ds_outputs['X']
    data_map = ds_outputs['data_map']
    is_valid = ds_outputs['is_valid']
    Y_ext = ds_outputs['Y_ext']
    V_ext = ds_outputs['V_ext'].copy()
    V = ds_outputs['V'].copy()

    m_map = {m['name']: m for m in src_models}
    rule_names = [r['name'] for r in a_rules]

    # 2.1 A-Rules
    print("\n\033[1m--- 2.1 A-Type Rules ---\033[0m\n")
    penalty_rows = []
    for col in tqdm(Y_ext.columns, desc="Applying A-Rules"):
        model_name, y_param = col.rsplit('__', 1)
        row = {'Model Name': model_name, 'Param': y_param}
        active_mask = Y_ext[col].notnull().values & is_valid
        N_active = np.sum(active_mask)
        if N_active < 30:
            for rule in a_rules:
                row[rule['name']] = "X"
            penalty_rows.append(row)
            continue
        model = m_map.get(model_name)
        if not model:
            continue
        model_context = model.get('inContext', []) + model.get('outContext', [])
        for rule in a_rules:
            if not set(rule['context']).issubset(set(model_context)) or Y_ext[col].isnull().all():
                row[rule['name']] = "X"
                continue
            new_scores = calculate_a_rule_validity(rule, active_mask, model, data_map, Y_ext, X)
            if np.all(np.isnan(new_scores)):
                row[rule['name']] = "X"
            else:
                mean_score = np.nanmean(new_scores)
                penalty = (1.0 - mean_score) * 100
                row[rule['name']] = round(penalty, 1)
            V_ext[col] = np.where(~np.isnan(V_ext[col].values), np.fmin(V_ext[col].values, new_scores), np.nan)
        penalty_rows.append(row)

    df_penalties = pd.DataFrame(penalty_rows)
    print("\n\033[1mA-Rule Physical Penalty Report (%):\033[0m")
    print(df_penalties.to_string(index=False))
    print("\nNote: A-Rules with many skipped models (X) should be converted into B-Rules.")

    # Plot validity distributions
    print("\nGenerating validity plots...")
    for out_param in output_params:
        p_name = out_param['name']
        Y_sub = Y_ext[[c for c in Y_ext.columns if c.endswith(f'__{p_name}')]]
        V_sub = V_ext[[c for c in V_ext.columns if c.endswith(f'__{p_name}')]]
        if not Y_sub.empty:
            plot_ensemble_validity(V_sub, Y_sub, out_param['desc'], output_dir='data')
        else:
            print(f"No data found for {out_param['desc']}")

    # 2.2 Initial Consensus
    print("\n\033[1m--- 2.2 Initial Consensus ---\033[0m\n")
    Y = pd.DataFrame(index=X.index)
    for param in output_params:
        p_name = param['name']
        p_cols = [col for col in Y_ext.columns if col.endswith(f'__{p_name}')]
        vals = Y_ext[p_cols].values
        wgts = V_ext[p_cols].values
        cons, conf = weighted_median(vals, wgts)
        Y[p_name] = cons
        V[p_name] = conf
    print(f"\033[32mInitial consensus calculated for {len(output_params)} parameters.\033[0m\n")
    total_samples = len(X)
    print("\033[1mDesign Space Occupancy:\033[0m")
    for col in Y.columns:
        count = Y[col].notna().sum()
        print(f"  Design space coverage for {col:<12} is {(count/total_samples*100):.2f}%")

    # 2.3 B-Rules
    print("\n\033[1m--- 2.3 B-Type Rules ---\033[0m\n")
    V_b = V_ext.copy()
    for rule in b_rules:
        V_b = calculate_b_rule_validity(X, Y, V_b, Y_ext, V_ext)
        print(f"B-Rule '{rule['name']}' applied.")

    # 2.4 Final Consensus
    print("\n\033[1m--- 2.4 Final Consensus ---\033[0m\n")
    for param in output_params:
        p_name = param['name']
        p_cols = [col for col in Y_ext.columns if col.endswith(f'__{p_name}')]
        vals = Y_ext[p_cols].values
        wgts = V_b[p_cols].values
        cons, conf = weighted_median(vals, wgts)
        Y[p_name] = cons
        V[p_name] = conf

    V_ext_aligned = check_model_alignment(Y_ext, Y, V_b)
    for param in output_params:
        p_name = param['name']
        p_cols = [col for col in Y_ext.columns if col.endswith(f'__{p_name}')]
        vals = Y_ext[p_cols].values
        wgts = V_ext_aligned[p_cols].values
        cons, conf = weighted_median(vals, wgts)
        Y[p_name] = cons
        V[p_name] = conf

    print("\nGenerating post-consensus validity plots...")
    for out_param in output_params:
        p_name = out_param['name']
        Y_sub = Y_ext[[c for c in Y_ext.columns if c.endswith(f'__{p_name}')]]
        V_sub = V_ext_aligned[[c for c in V_ext_aligned.columns if c.endswith(f'__{p_name}')]]
        if not Y_sub.empty:
            plot_ensemble_validity(V_sub, Y_sub, out_param['desc'] + " (Post-Consensus)", output_dir='data')
        else:
            print(f"No data found for {out_param['desc']} (Post-Consensus)")

    print(f"\033[32mFinal consensus calculated for {len(output_params)} parameters.\033[0m\n")
    print("\033[1mDesign Space Occupancy:\033[0m")
    total_samples = len(X)
    for col in Y.columns:
        count = Y[col].notna().sum()
        print(f"  Design space coverage for {col:<12} is {(count/total_samples*100):.2f}%")

    # Binary coverage matrix: ≥0.2 → covered
    covered = (V_ext_aligned >= 0.2).astype(int)
    n_models_covering = covered.sum(axis=1)  # total models covering each design point




    w_name, w_val = 25, 12
    bin_umbrella = f"{'BINARY COVERAGE':^{w_val * 5 + 4}}" 
    weight_umbrella = f"{'WEIGHTED COVERAGE':^{w_val * 2 + 1}}"
    umbrella_header = f"{'':<{w_name}}   {'':>{w_val}}   {bin_umbrella}   {weight_umbrella}"
    headers = (f"{'MODEL':<{w_name}} | {'CLAIMED':>{w_val}} |"
           f" {'SOLE':>{w_val}} {'CONTRIBUTING':>{w_val}} {'REDUNDANT':>{w_val}}"
           f" {'TOTAL':>{w_val}} {'REDUCTION':>{w_val}} |"
           f" {'COVERAGE':>{w_val}} {'REDUCTION':>{w_val}}")

    print(umbrella_header)
    print(headers)
    print("-" * len(headers))

    for col in Y_ext.columns.tolist():
        model_name = re.sub(r"([a-zA-Z]+)(\d{4}).*", lambda m: f"{m.group(1).capitalize()} ({m.group(2)})", col)

        # Claimed coverage (initial, before any rules)
        model_samples = Y_ext[col].count()
        claimed_pct = model_samples / total_samples

        # Weighted coverage (sum of final validity scores / total)
        weighted_pct = V_ext_aligned[col].sum() / total_samples

        # Binary coverage categories (V_ext_aligned >= 0.2)
        model_covered = covered[col].astype(bool)
        n_other = n_models_covering - model_covered
        pct_sole = (model_covered & (n_other == 0)).sum() / total_samples
        pct_contrib = (model_covered & (n_other >= 1) & (n_other < penalty_threshold-1)).sum() / total_samples
        pct_redund = (model_covered & (n_other >= penalty_threshold-1)).sum() / total_samples
        pct_binary_total = pct_sole + pct_contrib + pct_redund

        # Reductions (relative to originally claimed coverage)
        binary_reduction = (claimed_pct - pct_binary_total) / claimed_pct
        weighted_reduction = (claimed_pct - weighted_pct) / claimed_pct

        print(f"{model_name:<{w_name}} | {claimed_pct:>{w_val}.2%} |"
            f" {pct_sole:>{w_val}.2%} {pct_contrib:>{w_val}.2%} {pct_redund:>{w_val}.2%}"
            f" {pct_binary_total:>{w_val}.2%} {binary_reduction:>{w_val}.1%} |"
            f" {weighted_pct:>{w_val}.2%} {weighted_reduction:>{w_val}.1%}")

    print(f"* The average validity score across all data points, treating validity as a continuous spectrum to reflect the overall quality of the dataset.")

    return {
        'Y': Y,
        'V': V,
        'V_ext': V_ext_aligned,
    }


if __name__ == '__main__':
    print("This module should be run via run.py")