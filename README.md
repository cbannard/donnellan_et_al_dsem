# 

A two-step Bayesian pipeline for modelling infant behaviour time-series data and identifying predictors of language outcomes (CDI scores).

---

## Overview

The pipeline fits a multilevel Vector Autoregression (VAR) model to repeated-measures behaviour data, then uses Bayesian model averaging to identify which lag predictors (and their directions) best predict CDI vocabulary scores.

```
Step 1:   fitVAR.R                               → dsem_step1_*.rds
Step 1b:  granger_analysis.R          (optional) → granger_loocv_*.csv/.rds
Step 2a:  select_terms_for_individual_behaviours.R → dsem_selection_*_single_behaviour.rds
Step 2b:  final_term_weighting.R                 → dsem_selection_from_single_behaviour_sets*.csv/.rds
```

---

## Scripts

### `fitVAR.R` — Step 1: Multilevel VAR

Fits a Bayesian multilevel logistic VAR(p) model to all behaviours jointly, using a non-centred parameterisation over participants. All behaviours are binary, so each person has lag coefficients and intercepts only — there are no residual-variance (log-variance) parameters, since a Bernoulli outcome's variance is fixed by its probability. The script stops if any behaviour is not 0/1. Initialises HMC with Pathfinder variational inference for faster warmup.

**Usage**
```bash
Rscript fitVAR.R [job_id]
```
`job_id` defaults to `Sys.getenv("SLURM_JOB_ID")` or `"local"`.

**Key parameters** (edit at top of script)

| Parameter | Default | Description |
|---|---|---|
| `T_subsample` | 3165 | Max time points per participant |
| `p` | 5 | VAR lag order |
| `grainsize` | 1 | Parallelism grain for `partial_sum_lpmf` |
| `chains` | 1 | HMC chains |
| `iter_warmup` | 200 | Warmup iterations |
| `iter_sampling` | 400 | Post-warmup draws |
| `adapt_delta` | 0.95 | HMC target acceptance rate |
| `use_pathfinder` | TRUE | Use Pathfinder to initialise HMC |

**Input files**
- `data_full/*.csv` — per-participant behaviour time series (auto-detected by pattern)
- `CDI_data.csv` — CDI vocabulary scores (not used directly here; required in Steps 2a/2b)

**Output**
- `dsem_step1_<job_id>.rds` — list with:
  - `$fit` — raw cmdstanr fit object
  - `$draws` — posterior draws (`posterior` package format), including `z_person[i,q]`
  - `$prep` — data preparation object (see [RDS structure](#rds-structure) below)

**Dependencies**
- `cmdstanr` (≥ 0.5.3 for Pathfinder), CmdStan ≥ 2.33
- `posterior`, `dplyr`, `tidyr`, `readr`

---

### `granger_analysis.R` — Step 1b: Granger-causality screening (optional)

For each directed pair (x → y), computes the PSIS-LOO ELPD difference between the full VAR model and a reduced model with x's lag coefficients zeroed out for outcome y. A positive difference supports Granger causality. This step is independent of Steps 2a/2b and its output is not consumed downstream.

**Usage**
```bash
# Single file
Rscript granger_analysis.R dsem_step1_1.rds

# Folder of files (combined automatically)
Rscript granger_analysis.R MODELS-STEP1/ granger_loocv

# SLURM array job (one shard per task)
Rscript granger_analysis.R MODELS-STEP1/ granger_loocv array_task_id=1 array_task_count=9
```

**Input**
- Step 1 RDS file or folder of RDS files

**Output**
- `<output_prefix>_pairwise_elpd.csv` and `.rds` — one row per directed pair (x, y) with columns:
  `predictor_x`, `outcome_y`, `elpd_full`, `elpd_reduced_drop_xlags`, `elpd_diff_full_minus_reduced`, `se_diff`, `z_score`, `supports_x_granger_causing_y`, `strong_support_diff_gt_2se`, `mean_coef_x_to_y`, `prop_positive_x_to_y`, `coef_direction_x_to_y`

**Dependencies**
- `posterior`, `loo`

---

### `select_terms_for_individual_behaviours.R` — Step 2a: Within-behaviour block selection

For each outcome behaviour in turn, fits an error-in-variables Gaussian bridge model of CDI scores on VAR parameters (posterior means and SDs as uncertain predictors). Enumerates all 2^(V+1) − 1 subsets of predictor blocks for that behaviour (V lag blocks plus the mean), computes PSIS-LOO for each, and derives per-block inclusion weights via stacking and pseudoBMA+.

**Usage**
```bash
Rscript select_terms_for_individual_behaviours.R <behaviour_name> [step1_folder]
```
Runs one behaviour per call (`step1_folder` defaults to `MODELS-DSEM`); results are written to `step1_folder`. Run as a SLURM array job with one task per behaviour.

**Input**
- `dsem_step1_*.rds` or `dsem_result_*.rds` in the step1 folder
- `CDI_data.csv` — must be present in the **current working directory**

**Output** (one file per behaviour in output folder)
- `dsem_selection_<behaviour>_single_behaviour.rds` — list with:
  - `$behaviour` — behaviour name (character)
  - `$block_inclusion` — data frame: `block_idx`, `block_label`, `inclusion_weight`, `effect_direction` (or `NULL` if no combos fitted)
  - `$combo_summary`, `$stacking_weights`, `$pseudobma_weights`, `$baseline_elpd`, `$n_cc`

**Dependencies**
- `cmdstanr`, `posterior`, `loo`, `dplyr`

---

### `final_term_weighting.R` — Step 2b: Two-level stacking across behaviours

Combines the within-behaviour block inclusion weights from Step 2a with outer stacking weights across behaviour-set combinations to produce final per-term inclusion weights.

The final weight for a predictor term `q` is:

```
w_final(q) = outer_inclusion(behaviour(q)) × inner_inclusion(q | behaviour)
```

where `outer_inclusion` sums stacking weights over all behaviour-set combinations containing that behaviour, and `inner_inclusion` is the within-behaviour block weight from Step 2a.

**Usage**
```bash
Rscript final_term_weighting.R [selection_folder] [step1_folder] [output_prefix]
```
Defaults: `selection_folder = "DSEM-MODELS-STEP1-PERCULATE"`, `step1_folder = selection_folder`, `output_prefix = "dsem_selection_from_single_behaviour_sets"`.

**Input**
- `dsem_selection_*_single_behaviour.rds` in the selection folder (Step 2a output)
- `dsem_step1_*.rds` or `dsem_result_*.rds` in the step1 folder (Step 1 output)
- `CDI_data.csv` — must be present in the **current working directory**
- Optional: a checkpoint file (`<output_prefix>_checkpoint.rds`) to resume a partial run

**Output** (all saved to `selection_folder`)
- `<output_prefix>_stacking_coefs.csv/.rds` — per-term inclusion weights
- `<output_prefix>_stacking_weights.csv/.rds` — per-behaviour-combination stacking weights
- `<output_prefix>.csv/.rds` — combined summary

**Dependencies**
- `cmdstanr`, `posterior`, `loo`, `dplyr`

---

## RDS structure

All Step 1 RDS files have the same structure, consumed by all downstream scripts:

```
dsem_step1_*.rds
├── $fit         — cmdstanr CmdStanMCMC object
├── $draws       — posterior draws_matrix (includes z_person[i,q])
└── $prep
    ├── $stan_data       — list passed to Stan:
    │   ├── N, V, K, T_max, T_i
    │   ├── Y[N, T_max, V]       — outcome array
    │   ├── X_lag[N, T_max, K]   — lagged predictor array
    │   ├── obs_mask[N, T_max]   — 1 = observed time point
    │   ├── is_binary[V]         — 1 = binary outcome
    │   ├── prior_sd_A, prior_sd_intercept, prior_sd_sigma
    │   └── grainsize            — parallel sum grain size
    ├── $participants    — character vector, length N
    ├── $behaviours      — character vector, length V
    ├── $lag_col_names   — character vector, length K
    ├── $X_lag_center, $X_lag_scale  — standardisation constants
    ├── $V               — number of behaviours
    ├── $K               — V × p (total lag columns)
    └── $Q               — V×K + V (parameters per person)
```

The `z_person[i,q]` parameters encode each participant's individual VAR parameters in a non-centred form. `Q = V×K + V` decomposes as:
- `q = 1 … V×K`: lag coefficients A[y,x,lag] (column-major: outcome varies fastest)
- `q = V×K+1 … V×K+V`: intercepts μ[y]

---

## Workflow

```bash
# 1. Fit the VAR model (SLURM array or single run)
Rscript fitVAR.R $SLURM_JOB_ID

# 1b. (Optional) Granger causality screening
Rscript granger_analysis.R MODELS-STEP1/ granger_loocv \
  array_task_id=$SLURM_ARRAY_TASK_ID array_task_count=9

# 2a. Within-behaviour block selection (one SLURM task per behaviour, or all at once)
Rscript select_terms_for_individual_behaviours.R gaze MODELS-DSEM/   # repeat for each behaviour

# 2b. Two-level stacking across behaviours
Rscript final_term_weighting.R MODELS-DSEM/ MODELS-DSEM/ dsem_selection_from_single_behaviour_sets
```

`CDI_data.csv` must be in the directory from which Steps 2a and 2b are run.

---

## R package dependencies

| Package | Used by |
|---|---|
| `cmdstanr` | fitVAR, select_terms, final_term_weighting |
| `posterior` | fitVAR, granger_analysis, select_terms, final_term_weighting |
| `loo` | granger_analysis, select_terms, final_term_weighting |
| `dplyr` | fitVAR, select_terms, final_term_weighting |
| `tidyr` | fitVAR |
| `readr` | fitVAR |

CmdStan ≥ 2.33 is required for Pathfinder initialisation in `fitVAR.R`.
