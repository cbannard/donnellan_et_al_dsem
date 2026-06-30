# DSEMCODE

Scripts implementing a two-step Dynamic Structural Equation Model (DSEM) pipeline for infant behaviour time-series data, with optional Granger-style analysis.

---

## Pipeline overview

**Step 1 — `fitVAR.R`**
Fits a multilevel VAR (Vector Autoregression) model to infant behaviour time-series data using Stan via cmdstanr. Each participant gets person-level lag coefficients (A), intercepts (mu), and log-variances (logvar), estimated with a non-centred hierarchical prior. The posterior draws are saved to an RDS file. Multiple independent runs (one per SLURM job) can be combined downstream.

**Step 2a — `select_terms_for_individual_behaviours.R`**
For a single target behaviour, enumerates all subsets of predictor blocks (lag blocks for each of the 9 behaviours, plus mean and log-variance terms), fits each combination as a regression on CDI vocabulary outcome (log1p-transformed), and computes PSIS-LOO ELPD. Stacking weights and pseudo-BMA+ weights are used to derive per-block inclusion weights and effect directions. Run once per behaviour (9 runs total, one per behaviour).

**Step 2b — `final_term_weighting.R`**
Loads the per-behaviour selection results from Step 2a and forms all combinations of the 9 behaviour sets (511 combinations). Fits each combination on the shared complete-case CDI sample and computes outer stacking weights via `loo::stacking_weights`. Final term inclusion weights are the product of outer (behaviour-level) and inner (block-level) weights. Supports checkpoint/resume: progress is saved every 25 combinations and automatically resumed if the script is restarted.

**Optional — `granger_analysis.R`**
Granger-style causal scoring. For each directed pair x → y, compares the full-model PSIS-LOO ELPD against a reduced model with all x-lag coefficients zeroed. The ELPD difference quantifies the predictive contribution of x to y under the fitted VAR posterior. Designed to run as a SLURM array job, with each task handling a round-robin subset of outcome variables.

---

## Required R packages

- `cmdstanr` (with a working CmdStan installation)
- `posterior`
- `loo`
- `dplyr`
- `parallel` (base R, used by `granger_analysis.R`)

---

## Required data files

| File | Used by | Contents |
|------|---------|----------|
| `data_full/` | `fitVAR.R` | Per-participant CSV files named `<id>_11m.csv` and `<id>_12m.csv`, each with 9 behaviour columns |
| `CDI_data.csv` | `select_terms_for_individual_behaviours.R` `final_term_weighting.R` | Vocabulary outcome data; must contain columns `participant` and `us18_imp` |

---

## Running the pipeline

### Step 1: Fit the VAR model

```bash
Rscript fitVAR.R <jobid>
```

- `<jobid>` is an integer identifier; used to seed the sampler and name the output file.
- Output: `dsem_step1_<jobid>.rds` containing `$fit`, `$draws`, and `$prep`.

For a SLURM array job, submit one task per desired chain group and pass `$SLURM_ARRAY_TASK_ID` as the jobid. Example SLURM header:

```bash
#SBATCH --array=1-20
Rscript fitVAR.R $SLURM_ARRAY_TASK_ID
```

All resulting RDS files should be placed in a single folder for Steps 2a and 2b to load and combine.

---

### Step 2a: Per-behaviour predictor selection

Run once for each of the 9 behaviours. These can be submitted as independent jobs in parallel.

```bash
Rscript select_terms_for_individual_behaviours.R <behaviour> [input_folder]
```

- `<behaviour>`: one of `contingent`, `gaze`, `CV`, `non_CV`, `index_point`, `open_point`, `give`, `show`, `other_gesture`
- `[input_folder]`: directory containing the Step 1 RDS files (default: `MODELS-AGGREGATE-MODELSONLY`)

Output (written into `input_folder`):
- `dsem_selection_<behaviour>_single_behaviour.rds` — per-block inclusion weights and effect directions

Example SLURM array:

```bash
#SBATCH --array=1-9
behaviours=(contingent gaze CV non_CV index_point open_point give show other_gesture)
Rscript select_terms_for_individual_behaviours.R ${behaviours[$SLURM_ARRAY_TASK_ID-1]} MODELS-STEP1
```

---

### Step 2b: Two-level term weighting

Run after all 9 Step 2a jobs have completed.

```bash
Rscript final_term_weighting.R [selection_folder] [step1_folder] [output_prefix]
```

- `[selection_folder]`: folder containing the `dsem_selection_*_single_behaviour.rds` files and where output is written (default: `DSEM-MODELS-STEP1-PERCULATE`)
- `[step1_folder]`: folder containing Step 1 RDS files; can be the same as `selection_folder`
- `[output_prefix]`: prefix for output files (default: `dsem_selection_from_single_behaviour_sets`)

Outputs (6 files written to `selection_folder`):

| File | Contents |
|------|----------|
| `<prefix>_stacking_coefs.csv` / `.rds` | Per-term inclusion weights and effect directions |
| `<prefix>_stacking_weights.csv` / `.rds` | Per-combination stacking and pseudo-BMA+ weights |
| `<prefix>.rds` | Full results object including all intermediate structures |
| `<prefix>.csv` | Combination-level summary table |

If the script is interrupted, a checkpoint file (`<prefix>_checkpoint.rds`) is saved every 25 combinations. Re-running the same command will resume automatically from the last checkpoint.

---

### Optional: Granger causal scoring

```bash
Rscript granger_analysis.R <step1_rds_or_folder> [output_prefix] [loo_cores=1] [pair_cores=1] [chunk_size=2000]
```

For a SLURM array job with, for example, 9 tasks (one per outcome variable):

```bash
#SBATCH --array=1-9
Rscript granger_analysis.R MODELS-STEP1 granger_loocv array_task_id=$SLURM_ARRAY_TASK_ID array_task_count=$SLURM_ARRAY_TASK_COUNT
```

Outputs per task (written beside the input file or in the input folder):
- `<prefix>_task001of009_pairwise_elpd.rds`
- `<prefix>_task001of009_pairwise_elpd.csv`

After all tasks complete, combine shards:

```bash
Rscript dsem_var_granger_loocv_array_combine.R <output_dir> granger_loocv 9
```

This produces `granger_loocv_pairwise_elpd.rds` and `.csv` with all directed pairs ranked by ELPD advantage.

---

## Execution order summary

```
fitVAR.R  (×N jobs, one per chain group)
    ↓
select_terms_for_individual_behaviours.R  (×9 jobs, one per behaviour)
    ↓
final_term_weighting.R  (single job)

granger_analysis.R  (optional; ×V array tasks, one per outcome)
    ↓
dsem_var_granger_loocv_array_combine.R  (optional; single job to merge shards)
```

---

## Notes

- The bridge step (within Steps 2a and 2b) uses posterior means and SDs of the per-person VAR parameter vectors (`z_person`) to propagate Step 1 uncertainty into the Step 2 regression, rather than carrying full posterior draws forward.
- All CDI outcome values are log1p-transformed before regression.
- Lagged predictors are z-scored within Step 1 data preparation; the centering and scaling vectors are stored in `$prep$X_lag_center` and `$prep$X_lag_scale` in the Step 1 RDS.
