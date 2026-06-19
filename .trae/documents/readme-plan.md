# README Plan for Optimization Framework for Packed Beds

## Summary

Create a comprehensive `README.md` for the Optimization Framework repository. The README will serve as the primary landing page on GitHub, providing an overview of the framework, explaining its two core algorithms (cross-validation and gap detection), listing all implemented pressure drop models, and presenting the key results. The paper (PDF) and its appendix will be referenced but not duplicated — the README will be a concise, self-contained technical summary.

## Current State

The repository (`c:\Users\nikit\Documents\Optimization Framework`) contains:

- **`run.py`** — entry point orchestrating a 3-section pipeline: Dataset Generation → Cross-Validation → Gap Detection
- **`src/config.py`** — defines 6 input parameters (D, L, D_p, φ_s, ε, U_s_G_in), 1 output parameter (ΔP), 47 source models (all semi-empirical), 2 A-rules (Laminar/Turbulent Stress Tests), and 0 B-rules
- **`src/dataset_generation.py`** — generates design space via Sobol sequence (2^18 = 262,144 points), evaluates all source models with validity masking, produces initial validity matrices
- **`src/cross_validation.py`** — applies A-rules (physics-based Spearman correlation tests in local neighborhoods), computes weighted-median consensus, applies B-rules (energy balance residuals), produces final consensus and validity scores
- **`src/gap_detection.py`** — uses decision tree classifiers to identify blind-spot regions in the design space where model coverage is insufficient
- **`src/models_pressure_drop.py`** — 47 pressure drop correlations (1937–2026), each with validity ranges, DOI references, and standardized friction factor formulation
- **`src/models_misc.py`** — auxiliary models (gas properties, water properties, liquid holdup, void fraction, heat & mass transfer)
- **`data/`** — contains the draft paper PDF, LaTeX sources, figures, and CSV/PNG output files
- **`CONTRIBUTING.md`** — contributor guidelines for adding new models
- **`.github/ISSUE_TEMPLATE/new_model_data.yml`** — issue template for model submissions

## Proposed README Structure

### 1. Title & Badges
- Project title: "Optimization Framework for Packed Beds"
- Badges: Python version, license (MIT), maybe a DOI badge once published

### 2. Overview (~1 paragraph + diagram)
- **What**: A computational framework for systematic cross-validation and gap analysis of packed bed correlations (currently focused on single-phase pressure drop)
- **Why**: Dozens of ΔP correlations exist with partially overlapping validity ranges; the framework answers: which models are reliable, where do they agree, and where are the blind spots in the design space?
- **How**: 3-stage pipeline — generate a comprehensive design space, cross-validate models against physical laws and each other, detect gaps in coverage
- Reference the pipeline diagram from `data/cross-validation-routine.pdf` or `data/consensus-routine.pdf` (or a simplified ASCII/mermaid version)

### 3. Pipeline Description

#### 3.1 Dataset Generation
- Sobol sequence (2^18 ≈ 262k points) over 6 input parameters: D, L, D_p, φ_s, ε, U_s_G_in
- Each model is evaluated at every point; validity masks (from original publication ranges) produce NaN for out-of-range predictions
- Output: input matrix X, extended output matrix Y_ext (one column per model), initial validity matrix V_ext

#### 3.2 Cross-Validation
- **A-Rules (Physics-Based)**: Test whether each model respects known physical laws within local neighborhoods in the design space. Uses Spearman correlation to verify expected monotonic relationships (e.g., ΔP/L should increase with viscous stress in laminar regime). Currently 2 rules: Laminar Stress Test and Turbulent Stress Test.
- **Weighted Median Consensus**: Compute a consensus ΔP value at each design point as the weighted median of all valid model predictions, weighted by their validity scores.
- **B-Rules (Global Consistency)**: Check model outputs against global physical constraints (e.g., energy balance). Currently no B-rules configured.
- **Model Alignment**: Identify models whose predictions deviate significantly from the consensus, applying penalties or discarding outliers.
- **Final Consensus**: Re-compute consensus after applying all penalties and alignment checks.

#### 3.3 Gap Detection
- Convert final validity scores into a binary coverage map (V ≥ 0.2 → covered)
- Train a decision tree classifier to partition the design space into regions of high vs. low coverage
- Extract interpretable rules describing "blind spots" — regions of the design space where no model achieves sufficient validity
- Report coverage statistics: percentage of design space covered at various validity thresholds

### 4. Source Models (ΔP Correlations)
- Table listing all 47 pressure drop correlations with:
  - Model name (e.g., "Carman (1937)")
  - Year
  - Input parameters used (Re_p, L, ε, D, φ_s)
  - Reference DOI link
- Grouped or sorted by year for readability
- Note: all models are semi-empirical

### 5. Key Results
- **Design space coverage**: ~X% of the design space has at least one model, ~Y% has ≥5 overlapping models
- **Model validity distribution**: summary of the `validity_Pressure_drop.png` output — which models have high/medium/low physics validity
- **Consensus coverage**: what fraction of the design space is covered by the consensus after cross-validation
- **Blind spots**: key regions identified where no model is reliable (e.g., certain D/D_p ratios, Re ranges)
- Reference the paper for full analysis; include the key figures inline (e.g., `coverage-histogram.png`, `validity_Pressure_drop.png`, `gap-detection.png`)

### 6. Installation & Quick Start
- Requirements: Python 3.x, numpy, pandas, scipy, scikit-learn, bottleneck, tqdm, matplotlib
- Installation: `git clone` + `pip install -r requirements.txt` (create if missing)
- Quick start: `python run.py`

### 7. Repository Structure
- Directory tree showing `src/`, `data/`, `.github/`

### 8. Contributing
- Link to `CONTRIBUTING.md`
- Mention the issue template for new model submissions

### 9. License
- MIT License

### 10. Citation
- Citation block for the paper (once published)

## Implementation Steps

1. **Create `requirements.txt`** if missing — list dependencies: numpy, pandas, scipy, scikit-learn, bottleneck, tqdm, matplotlib
2. **Write `README.md`** following the structure above, using markdown formatting
3. **Verify** all referenced file paths and image links are correct

## Assumptions & Decisions

- The README will reference the paper PDF in `data/` but not reproduce its full content — the README is a technical summary for GitHub visitors
- Figure references will use relative paths to files in `data/` (e.g., `data/validity_Pressure_drop.png`)
- The model list will be generated from the actual `src_models` list in `config.py` and `__all__` in `models_pressure_drop.py` to ensure accuracy
- A `requirements.txt` will be created since one doesn't exist yet
- The README will use the "Erdim (2015)" style naming convention already established in the codebase (`pretty_label` function)

## Verification

- After writing README.md, review it renders correctly on GitHub (markdown preview)
- Verify all image paths resolve to actual files in `data/`
- Verify all DOI links are functional
- Ensure the model list matches the actual 47 models in `config.py`