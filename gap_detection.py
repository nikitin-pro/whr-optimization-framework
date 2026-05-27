import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeClassifier, _tree

from config import output_params


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
                    if vmin != float('-inf'):
                        rule_parts.append(f"{feat} > {vmin:.3f}")
                    if vmax != float('inf'):
                        rule_parts.append(f"{feat} <= {vmax:.3f}")
                blind_spots.append({
                    'Size': int(samples),
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

    return cv_outputs


if __name__ == '__main__':
    print("This module should be run via run.py")