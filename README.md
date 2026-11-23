# Study 1 (DSEM): temporal organisation of infant–caregiver interaction and later vocabulary

Code and data for Study 1 of *Identifying interaction sequences that drive
language development* (Phil. Trans. R. Soc. B, RSTB-2025-0407). Study 2 (the
CNN) has its own repository, linked from the paper.

Step 1 fits a Bayesian multilevel logistic vector autoregression (VAR(5)) to the
coded behaviour time series of 58 infant–caregiver dyads and tests which
behaviours predict which (Granger-style PSIS-LOO comparisons; Figure 1,
Table S1). Step 2 asks which features of each child's estimated dynamics predict
18-month expressive vocabulary: each predictor on its own (single-predictor
test; Figure 2), behaviour groups in a model search with pseudo-BMA+ weights,
and permutation tests against shuffled vocabulary (Table S2 and the text).

---

## Versions

The `main` branch reproduces the revised paper: the session-aware model
(`STEP1_SESSION=1`) and all of the permutation tests, with the settings given
under *Reproducing the paper*. The code of the original submission is on the
`original_submission` branch. With no options set, the scripts on `main`
reproduce the earlier analysis in which the two recording sessions were joined
end to end (commit `9b37810`); it is not reported in the revised paper.

---

## Data

- **`data_full.zip`**: the coded behaviour time series, one CSV per child and
  session, `data_full/p<child>_11m.csv` and `p<child>_12m.csv` (116 files,
  58 children × 2 sessions). Each file has a header row of time stamps
  (0, 0.2, 0.4, … s; 4,750 columns = 950 s at 5 Hz, shorter recordings padded
  with zeros) and one row per coded category, 0/1 for each 200 ms sample. The
  first nine rows are, in order, `contingent` (caregiver semantically
  contingent talk), `gaze` (infant gaze to the caregiver's face), `CV`,
  `non_CV` (vocalisations), `index_point`, `open_point`, `give`, `show` and
  `other_gesture`; the two remaining rows are not used. The coding scheme is
  described in the paper.
- **`CDI_data.csv`**: `participant`, `us18_imp`, the expressive vocabulary
  (Lincoln CDI) at 18 months, with missing values imputed from all available
  CDI measures by predictive mean matching (`mice`), as described in the paper.
  Step 2 models log(1 + `us18_imp`).

The raw video recordings are not shared (see the paper's Data Availability
Statement).

---

## Requirements

- A Slurm cluster. `run_dsem_pipeline.slurm` submits every stage as its own job
  with dependencies. Its settings (partition `himem`, module
  `apps/gcc/R/4.4.2`) are those of the University of Manchester CSF3 and can be
  changed at the top of the script (`PARTITION`, `R_MODULE`).
- R 4.4.2 with `cmdstanr`, `posterior`, `loo`, `dplyr` and `ggplot2`
  (optional: `ragg` for PNG figures, `bayesplot` for convergence plots), and
  CmdStan 2.40.0 (≥ 2.33 is needed for Pathfinder).
- Compute: Step 1 is 16 single-chain fits, each several hours to a day and
  about 1.7 GB of output. Granger and Step 2 take a few hours. Each group-level
  permutation analysis takes about 5.5 h of CPU per shuffle, the
  single-predictor test about 6 min per shuffle (one CPU per task).

---

## Reproducing the paper

Run from a login node in a checkout of `main`. Everything is written under
`RUN_DIR`; nothing is written to the repository. Do not edit or `git pull` the
scripts while jobs from the checkout are queued or running: jobs read the R
scripts when they start.

```bash
# Settings used for every result in the paper
export RUN_DIR=/scratch/$USER/dsem_session
export STEP1_SESSION=1                       # lags within each session + session intercepts
export STEP2_WEIGHTS=pseudobma STEP2B_BLOCKS=all
export FIG1_Z=z_child FIG2_WEIGHT=gain       # Figure 1: child-clustered z; Figure 2: single-predictor gain
export PERM_MAX_RUNNING=64                   # tasks of one array allowed to run at once

# 1. Step 1, convergence, Granger (Figure 1, Table S1), Step 2a/2b, figures, tables
bash run_dsem_pipeline.slurm

# 2. Single-predictor test: observed + 1,000 shuffles (Figure 2, q-values)
PERM_N=1000 bash run_dsem_pipeline.slurm single_block

# 3. Group-level permutation tests (Step 2b only). Set N to the number of
#    shuffles reported in the paper's Methods.
N=200
PERM_N=$N bash run_dsem_pipeline.slurm perm_groups                        # grouped by behaviour predicted
STEP2B_GROUPING=source PERM_N=$N bash run_dsem_pipeline.slurm perm_groups # grouped by preceding behaviour
STEP2B_EXCLUDE=contingent.contingent PERM_N=$N bash run_dsem_pipeline.slurm perm_groups
STEP2B_GROUPING=source STEP2B_EXCLUDE=contingent.contingent PERM_N=$N bash run_dsem_pipeline.slurm perm_groups
STEP2B_EXCLUDE=contingent.contingent,mean.contingent PERM_N=$N bash run_dsem_pipeline.slurm perm_groups
STEP2B_GROUPING=source STEP2B_EXCLUDE=contingent.contingent,mean.contingent PERM_N=$N bash run_dsem_pipeline.slurm perm_groups

# 4. Joint tests of contingent responses to vocalisations and to all behaviours
#    (1,000 shuffles; submit from the repository folder)
sbatch --array=0-10 -o $RUN_DIR/logs/block_sets_%a.out run_block_sets.slurm
#    when all 11 tasks have finished:
Rscript single_block_test.R summary $RUN_DIR sets=x set_tag=responses

# 5. Redraw Figure 2 and the tables once everything has finished
FIGURE=2 bash run_dsem_pipeline.slurm figures
bash run_dsem_pipeline.slurm tables
```

Stages 2–4 can be submitted together after stage 1 has finished; they wait for
the bridge file (`PERMUTATIONS/bridge.rds`, a per-child summary of the Step 1
posterior), which the first of them creates. A summary job runs after each
array. To re-run one by hand, use `Rscript permutation_null.R $RUN_DIR <prefix>`,
where `<prefix>` is `dsem_selection_from_single_behaviour_sets` followed by the
analysis's tag, e.g. `dsem_selection_from_single_behaviour_sets_src_excl_contingent_contingent`.
Shuffles can be added to any analysis later with the same settings and
`PERM_START=<first new> PERM_N=<last>`; finished shuffles are skipped.

**Seeds.** Step 1 chain *n* uses seed *n*. Shuffle *k* permutes the vocabulary
scores with seed `PERM_SEED_BASE + k` (default 20261112), the same in every
permutation analysis, and fits its models with seed `PERM_SEED_BASE + 1000k`.
Observed Step 2 runs draw a random seed unless `STEP2_SEED` is set; the seed is
printed in each log (`logs/step2a_*.out`, `logs/step2b.out`) and stored in the
output, so a run can be repeated exactly by passing `seed=<value>` to the
script.

### Where each result comes from (all under `RUN_DIR`)

| Paper | File |
|---|---|
| Convergence (Methods) | `CONVERGENCE/convergence_report.txt` |
| Figure 1, Step 1 results, Table S1 | `FIGURES/figure1_*`, `GRANGER/granger_loocv_pairwise_elpd.csv`, `GRANGER/granger_direction.csv`, `TABLES/table_S1_granger_pairs.csv`; Pareto-k counts in `logs/granger_*.out` |
| Figure 2 and the single-predictor results (gain, p, corrected p, q, effects) | `FIGURES/figure2_*`, `PERMUTATIONS/single_block_null.csv`, `PERMUTATIONS/single_block_report.txt` |
| Group and within-group inclusion probabilities, model ELPDs | `MODELS-STEP1/dsem_selection_from_single_behaviour_sets*_stacking_weights.csv`, `TABLES/table_S2_step2_blocks.csv` |
| Group-level p-values (main analyses and reduced groups) | `PERMUTATIONS/permutation_null[_<tag>]_report.txt`, tags `src`, `excl_contingent_contingent`, `src_excl_contingent_contingent`, `excl_contingent_contingent_mean_contingent`, `src_excl_contingent_contingent_mean_contingent` |
| Joint tests of contingent responses | `PERMUTATIONS/block_sets_responses_report.txt` |
| Table S2 | `TABLES/table_S2_step2_blocks.csv` (the `final_inclusion`, `p_inner`, `p_final` and `p_final_max` columns are not used for this model) |

---

## Scripts

| Script | Stage | What it does |
|---|---|---|
| `run_dsem_pipeline.slurm` | all | Submits the stages below as Slurm jobs with dependencies; `bash run_dsem_pipeline.slurm <stage>` re-runs one stage. Options are environment variables (below, and in the comments at the top of the script). |
| `fitVAR.R` | Step 1 | Bayesian multilevel logistic VAR(5) in Stan (cmdstanr, Pathfinder initialisation), one chain per task. Child-specific intercepts and lag effects (non-centred). Priors: N(0, 5) population intercepts with a half-N(0, 2) child SD; N(0, 1) lags and session shifts with half-N(0, 1) child SDs. Lagged predictors are 0/1, centred, not scaled. Every third sample is predicted (3,165 per child). |
| `check_convergence.R` | Step 1 | R-hat, bulk/tail ESS, divergences, tree depth and E-BFMI across chains. |
| `granger_analysis.R` | Step 1b | For each of the 72 directed pairs x → y, the PSIS-LOO ELPD of the full model for y minus that of the model with x's lags set to zero, with observation-level, child-clustered (used in the paper) and block (20 modelled samples ≈ 12 s) SEs. |
| `granger_direction.R` | Step 1c | Posterior probability that the population-level x → y effect, summed over lags, is positive. |
| `select_terms_for_individual_behaviours.R` | Step 2a | Within one behaviour group, all 1,023 combinations of its 10 predictor blocks (9 lag blocks and the mean), each an error-in-variables regression of log vocabulary on the children's posterior means with their posterior SDs; PSIS leave-one-child-out; within-group inclusion probabilities. |
| `final_term_weighting.R` | Step 2b | All 511 combinations of behaviour groups and the intercept-only model; pseudo-BMA+ weights; group inclusion probabilities. Options include `blocks=all`, `weights=pseudobma`, `grouping=source` (group blocks by the predicting behaviour), `exclude=x.y[,…]` or `mean.y` (leave blocks out of every model), `bridge=` and `permute_seed=`. Resumes from its checkpoint if interrupted. |
| `make_bridge.R` | permutations | Per-child posterior mean and SD of every Step 2 predictor, so that permutation runs need not load the Step 1 draws. |
| `permutation_null.R` | permutations | Compares an observed Step 2b run with its shuffled runs: group p-values (`p`, and `p_max` corrected for the nine groups), the intercept-only model's weight and the best model's ELPD. |
| `single_block_test.R` | single-predictor test | Each of the 90 blocks alone vs the intercept-only model (ELPD gain), for the observed and the shuffled vocabulary scores; `summary` gives p, the familywise p over 90 blocks, permutation q-values, and the number and sum of positive gains. `sets=name=x.y+x.y;…` tests sets of blocks as single models instead. |
| `run_block_sets.slurm` | joint tests | Array job for `single_block_test.R sets=…` (task 0 observed, tasks 1–10 shuffles 1–1,000). |
| `make_figures.R` | figures | Figure 1 (`z=z_child`: strong edges solid) and Figure 2 (`weight=gain`: every predictor with a positive single-predictor gain; width = gain, colour = sign of its effect). |
| `make_tables.R` | tables | Tables S1 and S2. |
| `monitor_progress.sh` | — | Progress of running Step 2 jobs. |

### Pipeline options (environment variables)

| Variable | Paper | Meaning |
|---|---|---|
| `RUN_DIR` | — | Output folder |
| `STEP1_SESSION` | `1` | Lags within each session; session-specific intercepts |
| `STEP1_ONSET`, `STEP1_LAG_STEP` | `0`, `1` | Other Step 1 variants (not reported) |
| `STEP2_WEIGHTS` | `pseudobma` | Model weights in Step 2 |
| `STEP2B_BLOCKS` | `all` | Step 2b uses every block of a group |
| `STEP2B_GROUPING` | `target`, `source` | Group blocks by the predicted or by the preceding behaviour |
| `STEP2B_EXCLUDE` | see above | Blocks left out of every model (several: comma-separated) |
| `PERM_N`, `PERM_START` | see above | Shuffles to run |
| `PERM_CHAINS`, `PERM_ITER` | `2`, `500` | MCMC per model in shuffled runs (observed runs use 4 chains) |
| `PERM_SEED_BASE` | `20261112` | Shuffle seeds |
| `GRANGER_BLOCK_ROWS` | `20` | Block length for the block SE |
| `FIG1_Z`, `FIG2_WEIGHT` | `z_child`, `gain` | Figure settings |

### Running without Slurm

Each stage is a single Rscript call; the stage blocks near the end of
`run_dsem_pipeline.slurm` give the exact arguments. In order:
`fitVAR.R <chain>` (16 chains, with `STEP1_SESSION=1` set) →
`check_convergence.R` → `granger_analysis.R` → `granger_direction.R` →
`select_terms_for_individual_behaviours.R <behaviour>` (9 behaviours) →
`final_term_weighting.R` → `make_bridge.R` → for each shuffle,
`single_block_test.R fit` and `final_term_weighting.R … permute_seed=…` →
`permutation_null.R`, `single_block_test.R summary` → `make_figures.R`,
`make_tables.R`. `select_terms_for_individual_behaviours.R`,
`final_term_weighting.R` and `single_block_test.R` read `CDI_data.csv` from
the working directory.

---

## Step 1 output structure

```
MODELS-STEP1/dsem_step1_<chain>.rds
├── $fit, $draws   cmdstanr fit and posterior draws (incl. z_person[i,q])
└── $prep          stan_data, participants, behaviours, lag_col_names,
                   V (9), K (V×5), Q (V×K + V), options (session, onset, lag_step)
```
In `z_person[i, q]`, `q = 1 … V·K` are the lag coefficients A[y, x, lag]
(outcome varies fastest) and `q = V·K+1 … V·K+V` the intercepts.
