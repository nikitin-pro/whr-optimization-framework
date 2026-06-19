import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeClassifier, _tree

from .config import input_params, output_params

# Display info for dependent input params (ratios relative to their base)
# 'transform': 'reciprocal' inverts the value (x → 1/x) and flips inequality operators
_param_display = {}
for p in input_params:
    if p['type'] == 'dependent' and p['base'] == 'D':
        if p['name'] == 'D_p':
            _param_display[p['name']] = {'label': 'D/D_p', 'transform': 'reciprocal'}
        elif p['name'] == 'L':
            _param_display[p['name']] = {'label': 'L/D'}


def detect_blind_spot_specs(X, V, target_col, threshold=0.2):
    y = (V[target_col].isna() | (V[target_col] < threshold)).astype(int)
    if y.nunique() <= 1:
        if y.iloc[0] == 1:
            print(f"  \033[31mcompletely blind\033[0m")
        else:
            print(f"  \033[32mperfectly covered\033[0m")
        return pd.DataFrame()

    dt = DecisionTreeClassifier(max_leaf_nodes=32, min_samples_leaf=0.01, class_weight='balanced')
    dt.fit(X, y)

    tree_ = dt.tree_
    feature_names = X.columns
    blind_spots = []

    def extract_rules(node, bounds):
        if tree_.feature[node] == _tree.TREE_UNDEFINED:
            samples = tree_.n_node_samples[node]
            values = tree_.value[node][0]
            blind_count = values[1] if len(values) > 1 else (values[0] if y.iloc[0] == 1 else 0)
            confidence = blind_count / sum(values)
            if confidence > 0.5:
                rule_parts = []
                for feat, (vmin, vmax) in bounds.items():
                    info = _param_display.get(feat)
                    label = info['label'] if info else feat

                    if info and info.get('transform') == 'reciprocal':
                        # Stored value is D_p/D; display as D/D_p = 1/(D_p/D)
                        # Invert bounds and flip operators:
                        #   D_p/D > vmin  →  D/D_p < 1/vmin
                        #   D_p/D <= vmax →  D/D_p >= 1/vmax
                        if vmax != float('inf'):
                            rule_parts.append(f"{label} >= {1.0/vmax:.2f}")
                        if vmin != float('-inf'):
                            rule_parts.append(f"{label} < {1.0/vmin:.2f}")
                    else:
                        if vmin != float('-inf'):
                            rule_parts.append(f"{label} > {vmin:.2f}")
                        if vmax != float('inf'):
                            rule_parts.append(f"{label} <= {vmax:.2f}")
                blind_spots.append({
                    'Size': f"{samples / len(X):.1%}",
                    'Confidence': f"{confidence:.1%}",
                    'Blind Spot Ranges': " AND ".join(rule_parts)
                })
            return
        name = feature_names[tree_.feature[node]]
        threshold_val = tree_.threshold[node]
        left_bounds = {k: list(v) for k, v in bounds.items()}
        left_bounds[name][1] = min(left_bounds[name][1], threshold_val)
        extract_rules(tree_.children_left[node], left_bounds)
        right_bounds = {k: list(v) for k, v in bounds.items()}
        right_bounds[name][0] = max(right_bounds[name][0], threshold_val)
        extract_rules(tree_.children_right[node], right_bounds)

    initial_bounds = {name: [float('-inf'), float('inf')] for name in feature_names}
    extract_rules(0, initial_bounds)

    return pd.DataFrame(blind_spots).sort_values(by='Size', ascending=False)


def run(cv_outputs):
    print("\033[1m" + "=" * 60)
    print("SECTION 3: Gap Detection")
    print("=" * 60 + "\033[0m")

    X = cv_outputs['X']
    V = cv_outputs['V']

    print(f"\nTotal design points: {len(X)}")

    for param in output_params:
        print(f"\n\033[1mBLIND SPOTS FOR \033[34m{param['desc']}\033[0m\033[1m:\033[0m")
        with pd.option_context('display.max_colwidth', None):
            blinds = detect_blind_spot_specs(X, V, param['name'])
            if not blinds.empty:
                print(blinds.head())
            else:
                print("  No significant blind spots detected.")

    # Design Space Coverage (validity bin distribution)
    print("\n\n\033[1mDESIGN SPACE COVERAGE:\033[0m")
    bins = np.arange(0, 1.1, 0.1)
    bin_labels = [f"{bins[i]:.1f}-{bins[i+1]:.1f}" for i in range(len(bins) - 1)]
    col_width = 11
    separator = "-" * (col_width * len(bin_labels) + 11)

    for param in output_params:
        p_name = param['name']
        values = V[p_name].dropna().values
        if len(values) == 0:
            continue

        counts, _ = np.histogram(values, bins=bins)
        percentages = (counts / len(X)) * 100

        print(f"\n\033[1m{param['desc']}:\033[0m")
        print(separator)
        header = "|" + "|".join(f" {label:<{col_width - 2}} " for label in bin_labels) + "|"
        print(header)
        print(separator)
        row = "|" + "|".join(f" {pct:<{col_width - 2}.1f} " for pct in percentages) + "|"
        print(row)
        print(separator)

        # Summary metrics
        avg_validity = np.nanmean(V[p_name].values)
        covered = (V[p_name].notna() & (V[p_name] >= 0.2)).sum() / len(X) * 100
        print(f"  Average validity: {avg_validity:.3f}  |  Design space with V ≥ 0.2: {covered:.1f}%")

    return cv_outputs


if __name__ == '__main__':
    print("This module should be run via run.py")