#!/bin/bash --login
# =============================================================================
# run_dsem_pipeline.slurm
#
# Runs the full Study 1 (DSEM) pipeline from the donnellan_et_al_dsem repo
# (Donnellan et al., Phil. Trans. R. Soc. B, RSTB-2025-0407):
#
#   Step 1   fitVAR.R                                 16 single-chain fits (1000 warmup + 1000 draws each)
#   Check    check_convergence.R                      R-hat, ESS, divergences, per-chain checks
#   Step 1b  granger_analysis.R                       Figure 1 / 72 directed pairs
#   Step 1c  granger_direction.R                      population-level direction of each pair
#   Step 2a  select_terms_for_individual_behaviours.R one task per behaviour
#   Step 2b  final_term_weighting.R                   Figure 2 inclusion weights
#   Figures  make_figures.R                           Figures 1 and 2 (network plots)
#   Tables   make_tables.R                            Supplementary tables S1 (Granger) and S2 (Step 2)
#
#   Permutation null for Step 2 (run separately, see below):
#   Bridge   make_bridge.R                            children's Step 1 posterior means/SDs, once
#   Perm 2a  select_terms_for_individual_behaviours.R PERM_N x 9 tasks, vocabulary shuffled
#   Perm 2b  final_term_weighting.R                   PERM_N tasks
#   Summary  permutation_null.R                       p-values against the permutation null
#
# HOW TO RUN
#   Put this file in the donnellan_et_al_dsem/ folder (next to fitVAR.R and
#   data_full.zip) on the cluster, then from a login node run:
#
#       bash run_dsem_pipeline.slurm
#
#   (use "bash", not "sbatch"). The script checks the setup, then submits every
#   stage as its own Slurm job, with dependencies so that each stage starts
#   only after the stages it needs have finished successfully:
#
#       step1 (array 1-16)  ──┬── convergence
#                             ├── granger (array) ── granger_combine ──┬── granger_direction
#                             │                                         │
#                             └── step2a  (array 0-8) ── step2b ───────┴── figures ── tables
#
#   To re-run a single stage (e.g. after fixing a failure):
#
#       bash run_dsem_pipeline.slurm step2a        # also: step1, convergence,
#                                               #   granger, granger_combine, granger_direction,
#                                               #   step2b, figures
#
#   Step 2 seeds: by default each Step 2 job draws a random seed and records
#   it (log + output .rds). To reproduce a run:  STEP2_SEED=12345 bash run_dsem_pipeline.slurm step2b
#
#   Step 2b keeps, for each behaviour, only blocks whose Step 2a inclusion
#   weight is >= 0.5. Change with e.g.  STEP2B_MIN_INNER=0.4 bash run_dsem_pipeline.slurm step2b
#   (or STEP2B_BLOCKS=best / all for the other block-selection modes).
#   Use pseudo-BMA+ instead of stacking weights at both levels with
#       STEP2_WEIGHTS=pseudobma STEP2B_MIN_INNER=0.6 bash run_dsem_pipeline.slurm step2b
#
#   Figures can be drawn one at a time:  FIGURE=1 bash run_dsem_pipeline.slurm figures
#   Figure 1 with the population-level positive rule (needs granger_direction):
#       FIGURE=1 FIG1_RULE=pop FIG1_PROB=0.95 bash run_dsem_pipeline.slurm figures
#
#   Step 1 skips chains whose dsem_step1_<n>.rds already exists, so
#   re-submitting step1 only reruns missing or failed chains.
#
# REVIEWER ROBUSTNESS ANALYSES
#   Step 1 options (see fitVAR.R; defaults reproduce the reported run). Use a
#   NEW RUN_DIR for each variant - the script refuses to mix options in one:
#     STEP1_SESSION=1   lags built within each session + session intercepts
#     STEP1_ONSET=1     predict behaviour onsets only (transition model)
#     STEP1_LAG_STEP=2  lags 400 ms apart (2 s window)
#   e.g. the recommended rerun:
#     RUN_DIR=/scratch/q55362cb/dsem_session STEP1_SESSION=1 \
#       STEP2_WEIGHTS=pseudobma STEP2B_BLOCKS=all bash run_dsem_pipeline.slurm
#     RUN_DIR=/scratch/q55362cb/dsem_onset STEP1_SESSION=1 STEP1_ONSET=1 \
#       bash run_dsem_pipeline.slurm step1      (then granger etc. as needed)
#
#   Granger SEs: granger_analysis.R now also reports child-clustered and block
#   SEs/z-scores (z_child, z_block). Re-score an existing Step 1 run without
#   refitting:   bash run_dsem_pipeline.slurm granger      (then granger_combine,
#   granger_direction, figures, tables). Block length: GRANGER_BLOCK_ROWS=20.
#   Figure 1 with strong (z > 2) edges solid and weak edges dashed, using the
#   child-clustered z:   FIGURE=1 FIG1_Z=z_child bash run_dsem_pipeline.slurm figures
#   (FIG1_SHOW=strong draws only the strong edges).
#
#   Figure 2 no longer uses the paper's arbitrary 1/3 cut-off: by default it
#   draws every block with final inclusion >= 0.01, line width = inclusion.
#   FIG2_THRESHOLD=perm draws only blocks above the permutation-calibrated
#   threshold (95th percentile of the largest final inclusion in the
#   permutation runs); FIG2_THRESHOLD=0.333 gives the old figure.
#   FIG2_WEIGHT=gain draws every block that on its own improves held-out
#   prediction of vocabulary (single-block test), line width = ELPD gain;
#   FIG2_WEIGHT=perm uses -log10(permutation p) instead. FIG2_MIN_GAIN=0.1
#   hides the weakest links.
#
#   Permutation null for Step 2 (after step2a + step2b have finished):
#     STEP2_WEIGHTS=pseudobma STEP2B_BLOCKS=all bash run_dsem_pipeline.slurm permutations
#   Use the SAME STEP2_WEIGHTS / STEP2B_BLOCKS / STEP2B_MIN_INNER as the observed
#   run. Submits bridge -> perm2a (PERM_N x 9 tasks) -> perm2b (PERM_N tasks)
#   -> perm_summary -> tables + Figure 2 (lines solid where the familywise
#   permutation p < .05, dashed otherwise). Settings: PERM_N (default 100), PERM_CHAINS (2)
#   and PERM_ITER (500) per model, PERM_MAX_RUNNING (100 tasks at once),
#   PERM_SEED_BASE. Each permutation costs roughly 9 x 2-3 h + 3-4 h of CPU;
#   finished permutations are skipped on resubmission, so it can be run in
#   batches (e.g. PERM_N=20 first, later PERM_N=100); to add permutations to an
#   existing set without queuing the finished ones, set PERM_START, e.g.
#   PERM_START=101 PERM_N=1000. Slurm limits array size (scontrol show config |
#   grep MaxArraySize); split big batches with PERM_START/PERM_N if needed.
#   Re-run only the summary
#   with:  bash run_dsem_pipeline.slurm perm_summary
#
#   Group-level permutations only (Step 2b without Step 2a; about 1/5 of the
#   cost). Gives the group and global p-values; block-level within-group and
#   final p-values keep using the runs that also have Step 2a. E.g. to go
#   from 100 to 1000 permutations:
#     STEP2_WEIGHTS=pseudobma STEP2B_BLOCKS=all PERM_START=101 PERM_N=1000 \
#       bash run_dsem_pipeline.slurm perm_groups
#
#   Variants of the group-level test (each with its own observed run, task 0,
#   and output files tagged in PERMUTATIONS/permutation_null_<tag>_*):
#     STEP2B_GROUPING=source   group blocks by the PREDICTING behaviour (how a
#                              behaviour organises what follows it)
#     STEP2B_EXCLUDE=contingent.contingent   leave a block out (here the
#                              contingent-talk self-lag)
#   e.g. STEP2_WEIGHTS=pseudobma STEP2B_BLOCKS=all STEP2B_EXCLUDE=contingent.contingent \
#          PERM_N=1000 bash run_dsem_pipeline.slurm perm_groups
#
#   Both group-level permutations and the single-block test in one go (the
#   bridge file is made first if it does not exist; everything else waits):
#     STEP2_WEIGHTS=pseudobma STEP2B_BLOCKS=all PERM_N=1000 bash run_dsem_pipeline.slurm group_tests
#
#   Single-block test (each predictor block alone vs intercept-only, with the
#   same shuffles; needs the bridge file, so run after/with permutations):
#     bash run_dsem_pipeline.slurm single_block
#   Submits PERM_N + 1 tasks (observed + shuffles, ~15 min each) -> summary
#   -> tables -> Figure 2 (solid = familywise single-block p < .05).
#   Results: PERMUTATIONS/single_block_null.csv, single_block_report.txt
#
# OUTPUT: everything the pipeline writes (results, logs, temp files and the
# unzipped data) goes under $RUN_DIR, default /scratch/q55362cb/dsem_results/.
# Nothing is written into the repo folder. Override with e.g.
#   RUN_DIR=/scratch/q55362cb/run2 bash run_dsem_pipeline.slurm
#   MODELS-STEP1/dsem_step1_<n>.rds                     Step 1 posterior draws
#   GRANGER/granger_loocv_pairwise_elpd.csv             Step 1b, all 72 pairs
#   MODELS-STEP1/dsem_selection_<b>_single_behaviour.rds  Step 2a
#   MODELS-STEP1/dsem_selection_from_single_behaviour_sets_stacking_coefs.csv
#                                                       Step 2b final weights
#   CONVERGENCE/convergence_report.txt                  Step 1 convergence summary (+ CSVs, plots)
#   GRANGER/granger_direction.csv                       pairs + population-level P(effect > 0)
#   FIGURES/figure1_granger_network.{pdf,png}            Figure 1 (+ figure1_edges.csv)
#   FIGURES/figure2_vocabulary_network.{pdf,png}         Figure 2 (+ figure2_terms.csv)
#   TABLES/table_S1_granger_pairs.csv, table_S2_step2_blocks.csv   supplementary tables
#   PERMUTATIONS/perm_<k>/                              permutation runs
#   PERMUTATIONS/permutation_null_*.csv, _report.txt    permutation p-values
#   logs/                                               one log per task
#
# Cluster settings below match the original runs (CSF3-style "himem"
# partition, module apps/gcc/R/4.4.2, 7-day limit). CmdStan (>= 2.33) and the
# R packages cmdstanr, posterior, loo and dplyr must already be
# installed for this R module.
# =============================================================================

# ---------------------------- configuration ---------------------------------
PARTITION="${PARTITION:-himem}"
R_MODULE="${R_MODULE:-apps/gcc/R/4.4.2}"

N_CHAINS="${N_CHAINS:-16}"           # Step 1: one chain per array task (seed = task id)
GRANGER_SHARDS="${GRANGER_SHARDS:-9}" # Step 1b: 9 = one outcome behaviour per task;
                                      #   1 = everything in one job (as originally run)

TIME_LONG="7-0"                      # wall time for model-fitting stages
MEM_STEP1="${MEM_STEP1:-}"           # blank = partition default (as originally run)
MEM_GRANGER="${MEM_GRANGER:-512G}"   # original Granger job used 512G
MEM_STEP2="${MEM_STEP2:-256G}"       # steps 2a/2b load all step-1 files at once
MEM_CONV="${MEM_CONV:-128G}"         # convergence check holds all chains' draws

STEP2B_BLOCKS="${STEP2B_BLOCKS:-min}"    # step2b: which Step 2a blocks enter (min | best | all)
STEP2B_MIN_INNER="${STEP2B_MIN_INNER:-0.5}"  # step2b: inclusion-weight cut-off for blocks=min
STEP2_WEIGHTS="${STEP2_WEIGHTS:-stacking}"     # step2b: stacking or pseudobma model weights
FIGURE="${FIGURE:-both}"             # figures stage: 1, 2 or both
FIG1_RULE="${FIG1_RULE:-mean}"       # Figure 1 positive rule: mean or pop
FIG1_PROB="${FIG1_PROB:-0.95}"       # probability cut-off when FIG1_RULE=pop
FIG1_Z="${FIG1_Z:-z_score}"          # Figure 1 z: z_score | z_child | z_block
FIG2_WEIGHT="${FIG2_WEIGHT:-inclusion}" # Figure 2 links: inclusion (two-level) | gain | perm
                                      #   (gain/perm: single-block evidence; needs single_block)
FIG2_MIN_GAIN="${FIG2_MIN_GAIN:-0}"  # Figure 2 with gain/perm: minimum ELPD gain to draw a link
FIG2_THRESHOLD="${FIG2_THRESHOLD:-0.01}" # Figure 2 blocks drawn: 0.01 = all with non-negligible
                                      #   weight; perm = permutation-calibrated threshold;
                                      #   or a number (the original paper used 0.333)
FIG1_SHOW="${FIG1_SHOW:-all}"        # Figure 1: all (weak edges dashed) | strong
STEP1_SESSION="${STEP1_SESSION:-0}"  # Step 1 options (fitVAR.R); defaults = reported run
STEP1_ONSET="${STEP1_ONSET:-0}"
STEP1_LAG_STEP="${STEP1_LAG_STEP:-1}"
GRANGER_BLOCK_ROWS="${GRANGER_BLOCK_ROWS:-20}"  # rows per block for the block SE
PERM_N="${PERM_N:-100}"              # permutation null: run permutations PERM_START..PERM_N
PERM_START="${PERM_START:-1}"        #   (e.g. PERM_START=101 PERM_N=1000 adds 900 to an existing 100)
STEP2B_GROUPING="${STEP2B_GROUPING:-target}" # perm_groups variants: target | source grouping
STEP2B_EXCLUDE="${STEP2B_EXCLUDE:-}"         #   blocks left out, e.g. contingent.contingent
PERM_CHAINS="${PERM_CHAINS:-2}"      #   MCMC chains per Step 2 model in permutation runs
PERM_ITER="${PERM_ITER:-500}"        #   warm-up (= sampling) iterations per chain
PERM_MAX_RUNNING="${PERM_MAX_RUNNING:-100}"  # max permutation tasks running at once
PERM_SEED_BASE="${PERM_SEED_BASE:-20261112}" # permutation k uses seed PERM_SEED_BASE + k
MEM_BRIDGE="${MEM_BRIDGE:-64G}"      # bridge job: reads one Step 1 chain at a time
MEM_PERM="${MEM_PERM:-16G}"          # permutation tasks read the small bridge file only
STEP2_SEED="${STEP2_SEED:-}"         # blank = each Step 2 job draws a random seed
                                      #   (printed in its log and saved in its .rds);
                                      #   set a number to reproduce an earlier run

BEHAVIOURS=(contingent gaze CV non_CV index_point open_point give show other_gesture)
# Variant group-level analyses write under their own Step 2b prefix.
VARIANT_TAG=""
[[ "${STEP2B_GROUPING:-target}" == source ]] && VARIANT_TAG="src"
[[ -n "${STEP2B_EXCLUDE:-}" ]] && VARIANT_TAG="${VARIANT_TAG:+${VARIANT_TAG}_}excl_${STEP2B_EXCLUDE//[.,]/_}"
GRANGER_PREFIX="granger_loocv"
STEP2B_PREFIX="dsem_selection_from_single_behaviour_sets"
VARIANT_PREFIX="${STEP2B_PREFIX}${VARIANT_TAG:+_$VARIANT_TAG}"
# -----------------------------------------------------------------------------

set -eo pipefail

STAGE="${STAGE:-${1:-submit}}"

# Paths. When submitting, REPO_DIR defaults to the folder holding this script;
# compute stages receive REPO_DIR/RUN_DIR through --export.
if [[ -z "${REPO_DIR:-}" ]]; then
  REPO_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
fi
RUN_DIR="${RUN_DIR:-/scratch/q55362cb/dsem_results}"
DATA_DIR="$RUN_DIR/data_full"   # data_full.zip is unzipped here, not in the repo
STEP1_DIR="$RUN_DIR/MODELS-STEP1"
GRANGER_DIR="$RUN_DIR/GRANGER"
PERM_DIR="$RUN_DIR/PERMUTATIONS"
BRIDGE_FILE="$PERM_DIR/bridge.rds"
LOG_DIR="$RUN_DIR/logs"

load_r() {
  module purge
  module load "$R_MODULE"
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
}

# Private temp dir per task: cmdstanr compiles models and writes CmdStan CSVs
# (~1 GB per Step 1 chain) to tempdir(), which can fill a node's /tmp.
use_private_tmp() {
  export TMPDIR="$RUN_DIR/tmp/${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
  mkdir -p "$TMPDIR"
  trap 'rm -rf "$TMPDIR"' EXIT
}

banner() {
  echo "=== $1 | job ${SLURM_JOB_ID:-n/a} task ${SLURM_ARRAY_TASK_ID:-n/a} | $(hostname) | $(date) ==="
}

# Builds one sbatch call for a stage. Extra sbatch flags come first.
submit_stage() {
  local stage="$1"; shift
  sbatch --parsable -p "$PARTITION" -n 1 \
    --export=ALL,STAGE="$stage",REPO_DIR="$REPO_DIR",RUN_DIR="$RUN_DIR",STEP2_SEED="$STEP2_SEED",FIGURE="$FIGURE",FIG1_RULE="$FIG1_RULE",FIG1_PROB="$FIG1_PROB",FIG1_Z="$FIG1_Z",FIG2_THRESHOLD="$FIG2_THRESHOLD",FIG2_WEIGHT="$FIG2_WEIGHT",FIG2_MIN_GAIN="$FIG2_MIN_GAIN",FIG1_SHOW="$FIG1_SHOW",STEP2B_BLOCKS="$STEP2B_BLOCKS",STEP2B_MIN_INNER="$STEP2B_MIN_INNER",STEP2_WEIGHTS="$STEP2_WEIGHTS",STEP1_SESSION="$STEP1_SESSION",STEP1_ONSET="$STEP1_ONSET",STEP1_LAG_STEP="$STEP1_LAG_STEP",GRANGER_BLOCK_ROWS="$GRANGER_BLOCK_ROWS",PERM_N="$PERM_N",PERM_START="$PERM_START",STEP2B_GROUPING="$STEP2B_GROUPING",STEP2B_EXCLUDE="$STEP2B_EXCLUDE",PERM_CHAINS="$PERM_CHAINS",PERM_ITER="$PERM_ITER",PERM_SEED_BASE="$PERM_SEED_BASE" \
    "$@" "$SELF"
}

mem_flag() { [[ -n "$1" ]] && echo "--mem=$1"; }

if [[ -n "${SLURM_JOB_ID:-}" && "$STAGE" == submit ]]; then
  echo "Run this file with 'bash', not 'sbatch' (see header)." >&2; exit 1
fi

# ----------------------------- submission ------------------------------------
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  SELF="$(readlink -f "${BASH_SOURCE[0]}")"
  WHAT="$STAGE"   # "submit" = whole pipeline; otherwise a single stage

  echo "Repo:    $REPO_DIR"
  echo "Run dir: $RUN_DIR"

  # --- preflight checks ---
  for f in fitVAR.R granger_analysis.R select_terms_for_individual_behaviours.R \
           final_term_weighting.R check_convergence.R make_figures.R granger_direction.R \
           make_tables.R make_bridge.R permutation_null.R single_block_test.R CDI_data.csv; do
    [[ -f "$REPO_DIR/$f" ]] || { echo "Missing $REPO_DIR/$f" >&2; exit 1; }
  done
  mkdir -p "$RUN_DIR" 2>/dev/null && [[ -w "$RUN_DIR" ]] \
    || { echo "Cannot create or write to $RUN_DIR" >&2; exit 1; }
  if [[ ! -d "$DATA_DIR" ]]; then
    [[ -f "$REPO_DIR/data_full.zip" ]] || { echo "Missing $REPO_DIR/data_full.zip" >&2; exit 1; }
    echo "Unzipping data_full.zip into $RUN_DIR ..."
    unzip -q "$REPO_DIR/data_full.zip" -x '__MACOSX/*' -d "$RUN_DIR"
  fi
  n_csv=$(find "$DATA_DIR" -name 'p*_1[12]m.csv' | wc -l)
  echo "Found $n_csv participant CSV files in data_full/ (expect 116 = 58 children x 2 sessions)"

  if command -v module >/dev/null 2>&1; then
    load_r >/dev/null 2>&1 || true
  fi
  if command -v Rscript >/dev/null 2>&1; then
    Rscript -e 'for (p in c("cmdstanr","posterior","loo","dplyr")) if (!requireNamespace(p, quietly=TRUE)) stop("R package missing: ", p); v <- cmdstanr::cmdstan_version(); cat("CmdStan", v, "\n"); if (package_version(v) < "2.33") stop("CmdStan >= 2.33 needed for Pathfinder")' \
      || { echo "Preflight R check failed (see above)." >&2; exit 1; }
  else
    echo "WARNING: Rscript not found on this node; skipping package check."
  fi

  mkdir -p "$STEP1_DIR" "$GRANGER_DIR" "$LOG_DIR" "$RUN_DIR/tmp"
  ln -sfn "$DATA_DIR" "$STEP1_DIR/data_full"   # fitVAR.R reads ./data_full

  # Step 1 options are recorded per run folder; refuse to mix them.
  opts="session=$STEP1_SESSION onset=$STEP1_ONSET lag_step=$STEP1_LAG_STEP"
  opts_file="$RUN_DIR/step1_options.txt"
  if [[ "$WHAT" == submit || "$WHAT" == step1 ]]; then
    if [[ -f "$opts_file" ]] && [[ "$(cat "$opts_file")" != "$opts" ]]; then
      echo "This RUN_DIR was used with Step 1 options '$(cat "$opts_file")', not '$opts'." >&2
      echo "Use a new RUN_DIR for a different Step 1 variant." >&2; exit 1
    fi
    if [[ ! -f "$opts_file" ]] && compgen -G "$STEP1_DIR/dsem_step1_*.rds" >/dev/null && [[ "$opts" != "session=0 onset=0 lag_step=1" ]]; then
      echo "$STEP1_DIR already holds Step 1 fits made without options; use a new RUN_DIR for '$opts'." >&2; exit 1
    fi
    echo "$opts" > "$opts_file"
    echo "Step 1 options: $opts"
  fi

  n_beh=${#BEHAVIOURS[@]}
  dep() { [[ "$WHAT" == submit && -n "$1" ]] && echo "--dependency=afterok:$1"; }

  J1=""; JG=""; J2A=""
  if [[ "$WHAT" == submit || "$WHAT" == step1 ]]; then
    J1=$(submit_stage step1 -J dsem_step1 -t "$TIME_LONG" $(mem_flag "$MEM_STEP1") \
           --array=1-"$N_CHAINS" -o "$LOG_DIR/step1_%a.out")
    echo "Step 1   (fitVAR, $N_CHAINS chains):           job $J1"
  fi
  if [[ "$WHAT" == submit || "$WHAT" == convergence ]]; then
    JC=$(submit_stage convergence -J dsem_convergence -t 1-0 $(mem_flag "$MEM_CONV") \
           $(dep "$J1") -o "$LOG_DIR/convergence.out")
    echo "Check    (Step 1 convergence):              job $JC"
  fi
  if [[ "$WHAT" == submit || "$WHAT" == granger ]]; then
    JG=$(submit_stage granger -J dsem_granger -t "$TIME_LONG" $(mem_flag "$MEM_GRANGER") \
           $(dep "$J1") --array=1-"$GRANGER_SHARDS" -o "$LOG_DIR/granger_%a.out")
    echo "Step 1b  (Granger, $GRANGER_SHARDS shards):          job $JG"
  fi
  if [[ "$WHAT" == submit || "$WHAT" == granger_combine ]]; then
    JGC=$(submit_stage granger_combine -J dsem_granger_combine -t 0-2 --mem=8G \
            $(dep "$JG") -o "$LOG_DIR/granger_combine.out")
    echo "Step 1b  (combine Granger shards):          job $JGC"
  fi
  if [[ "$WHAT" == submit || "$WHAT" == granger_direction ]]; then
    JGD=$(submit_stage granger_direction -J dsem_granger_direction -t 0-6 --mem=64G \
            $(dep "$JGC") -o "$LOG_DIR/granger_direction.out")
    echo "Step 1c  (population-level direction):      job $JGD"
  fi
  if [[ "$WHAT" == submit || "$WHAT" == step2a ]]; then
    J2A=$(submit_stage step2a -J dsem_step2a -t "$TIME_LONG" $(mem_flag "$MEM_STEP2") \
            $(dep "$J1") --array=0-$((n_beh - 1)) -o "$LOG_DIR/step2a_%a.out")
    echo "Step 2a  (per-behaviour selection, $n_beh tasks): job $J2A"
  fi
  if [[ "$WHAT" == submit || "$WHAT" == step2b ]]; then
    J2B=$(submit_stage step2b -J dsem_step2b -t "$TIME_LONG" $(mem_flag "$MEM_STEP2") \
            $(dep "$J2A") -o "$LOG_DIR/step2b.out")
    echo "Step 2b  (final term weighting):            job $J2B"
  fi
  if [[ "$WHAT" == submit || "$WHAT" == figures ]]; then
    fig_dep=""
    [[ "$WHAT" == submit ]] && fig_dep="--dependency=afterok:${JGC}:${J2B}"
    JF=$(submit_stage figures -J dsem_figures -t 0-1 --mem=8G $fig_dep -o "$LOG_DIR/figures.out")
    echo "Figures  (Figures 1 and 2):                 job $JF"
  fi
  if [[ "$WHAT" == submit || "$WHAT" == tables ]]; then
    tab_dep=""
    [[ "$WHAT" == submit ]] && tab_dep="--dependency=afterany:${JGD}:${J2B}"
    JT=$(submit_stage tables -J dsem_tables -t 0-1 --mem=8G $tab_dep -o "$LOG_DIR/tables.out")
    echo "Tables   (supplementary tables S1, S2):     job $JT"
  fi

  # ---- bridge file for the permutation stages: made first if missing ----
  # (bridge_dep echoes a --dependency flag on the bridge job, or nothing if
  # bridge.rds already exists; the bridge job is submitted at most once.)
  JB_AUTO=""
  ensure_bridge() {
    [[ -s "$BRIDGE_FILE" ]] && return 0
    mkdir -p "$PERM_DIR"
    if [[ -z "$JB_AUTO" ]]; then
      JB_AUTO=$(submit_stage bridge -J dsem_bridge -t 0-6 $(mem_flag "$MEM_BRIDGE") -o "$LOG_DIR/bridge.out")
      echo "Bridge   (Step 1 summaries; others wait for it): job $JB_AUTO"
    fi
  }
  bridge_dep() { [[ -n "$JB_AUTO" ]] && echo "--dependency=afterok:$JB_AUTO"; }
  if [[ "$WHAT" == group_tests || "$WHAT" == single_block || "$WHAT" == perm_groups ]]; then
    ensure_bridge
  fi

  # ---- single-block permutation test (needs the bridge file) ----
  if [[ "$WHAT" == single_block || "$WHAT" == sb_summary || "$WHAT" == group_tests ]]; then
    JSB=""
    if [[ "$WHAT" == single_block || "$WHAT" == group_tests ]]; then
      JSB=$(submit_stage sb_fit -J dsem_single_block -t 1-0 --mem=8G $(bridge_dep) \
              --array=$(( PERM_START == 1 ? 0 : PERM_START ))-"$PERM_N"%"$PERM_MAX_RUNNING" -o "$LOG_DIR/single_block_%a.out")
      echo "Single-block test (shuffles $PERM_START-$PERM_N$( (( PERM_START == 1 )) && echo ' + observed')): job $JSB"
    fi
    sdep=""; [[ -n "$JSB" ]] && sdep="--dependency=afterany:$JSB"
    JSS=$(submit_stage sb_summary -J dsem_sb_summary -t 0-1 --mem=8G $sdep -o "$LOG_DIR/single_block_summary.out")
    echo "Single-block summary:                       job $JSS"
    JT=$(submit_stage tables -J dsem_tables -t 0-1 --mem=8G --dependency=afterany:$JSS -o "$LOG_DIR/tables.out")
    echo "Tables:                                     job $JT"
    JF2=$(FIGURE=2 submit_stage figures -J dsem_figure2 -t 0-1 --mem=8G --dependency=afterok:$JSS \
            -o "$LOG_DIR/figure2_perm.out")
    echo "Figure 2:                                   job $JF2"
  fi

  # ---- group-level permutations: Step 2b only (no Step 2a) ----
  if [[ "$WHAT" == perm_groups || "$WHAT" == group_tests ]]; then
    [[ "$STEP2B_BLOCKS" == all ]] || { echo "perm_groups needs STEP2B_BLOCKS=all" >&2; exit 1; }
    mkdir -p "$PERM_DIR"
    first=$PERM_START
    # A variant (STEP2B_GROUPING=source or STEP2B_EXCLUDE) also needs its observed run (task 0).
    [[ -n "$VARIANT_TAG" && "$PERM_START" == 1 ]] && first=0
    JPG=$(submit_stage perm2b_only -J dsem_perm_groups -t "$TIME_LONG" --mem="$MEM_PERM" $(bridge_dep) \
            --array="$first"-"$PERM_N"%"$PERM_MAX_RUNNING" -o "$LOG_DIR/perm2b_only${VARIANT_TAG:+_$VARIANT_TAG}_%a.out")
    echo "Perm 2b only${VARIANT_TAG:+ [$VARIANT_TAG]} (runs $first-$PERM_N${VARIANT_TAG:+; 0 = observed}): job $JPG"
    JPS=$(submit_stage perm_summary -J dsem_perm_summary -t 0-2 --mem=8G --dependency=afterany:$JPG \
            -o "$LOG_DIR/perm_summary${VARIANT_TAG:+_$VARIANT_TAG}.out")
    echo "Summary  (permutation p-values):            job $JPS"
    if [[ -z "$VARIANT_TAG" ]]; then
      JT=$(submit_stage tables -J dsem_tables -t 0-1 --mem=8G --dependency=afterany:$JPS -o "$LOG_DIR/tables.out")
      echo "Tables:                                     job $JT"
    fi
  fi

  # ---- permutation null for Step 2 (not part of the default submission) ----
  if [[ "$WHAT" == permutations || "$WHAT" == bridge || "$WHAT" == perm2a || "$WHAT" == perm2b \
        || "$WHAT" == perm_summary ]]; then
    mkdir -p "$PERM_DIR"
    pdep() { [[ "$WHAT" == permutations && -n "$1" ]] && echo "--dependency=afterany:$1"; }
    JB=""; JP2A=""; JP2B=""
    if [[ "$WHAT" == permutations || "$WHAT" == bridge ]]; then
      if [[ "$WHAT" == bridge || ! -s "$BRIDGE_FILE" ]]; then
        JB=$(submit_stage bridge -J dsem_bridge -t 0-6 $(mem_flag "$MEM_BRIDGE") -o "$LOG_DIR/bridge.out")
        echo "Bridge   (Step 1 summaries for Step 2):     job $JB"
      fi
    fi
    if [[ "$WHAT" == permutations || "$WHAT" == perm2a ]]; then
      n_tasks=$(( (PERM_N - PERM_START + 1) * n_beh ))
      (( n_tasks > 0 )) || { echo "PERM_N must be >= PERM_START" >&2; exit 1; }
      bdep=""; [[ -n "$JB" ]] && bdep="--dependency=afterok:$JB"
      JP2A=$(submit_stage perm2a -J dsem_perm2a -t "$TIME_LONG" --mem="$MEM_PERM" $bdep \
               --array=0-$((n_tasks - 1))%"$PERM_MAX_RUNNING" -o "$LOG_DIR/perm2a_%a.out")
      echo "Perm 2a  (permutations $PERM_START-$PERM_N x $n_beh behaviours):  job $JP2A"
    fi
    if [[ "$WHAT" == permutations || "$WHAT" == perm2b ]]; then
      JP2B=$(submit_stage perm2b -J dsem_perm2b -t "$TIME_LONG" --mem="$MEM_PERM" $(pdep "$JP2A") \
               --array="$PERM_START"-"$PERM_N"%"$PERM_MAX_RUNNING" -o "$LOG_DIR/perm2b_%a.out")
      echo "Perm 2b  (permutations $PERM_START-$PERM_N):            job $JP2B"
    fi
    if [[ "$WHAT" == permutations || "$WHAT" == perm_summary ]]; then
      JPS=$(submit_stage perm_summary -J dsem_perm_summary -t 0-2 --mem=8G $(pdep "$JP2B") \
              -o "$LOG_DIR/perm_summary.out")
      echo "Summary  (permutation p-values):            job $JPS"
      if [[ "$WHAT" == permutations ]]; then
        JT=$(submit_stage tables -J dsem_tables -t 0-1 --mem=8G --dependency=afterany:$JPS -o "$LOG_DIR/tables.out")
        echo "Tables   (with permutation p-values):       job $JT"
        JF2=$(FIGURE=2 submit_stage figures -J dsem_figure2 -t 0-1 --mem=8G --dependency=afterok:$JPS \
                -o "$LOG_DIR/figure2_perm.out")
        echo "Figure 2 (with permutation significance):   job $JF2"
      fi
    fi
  fi
  echo "Monitor with: squeue -u \$USER    Logs: $LOG_DIR"
  exit 0
fi

# ------------------------------ compute stages -------------------------------
case "$STAGE" in

step1)
  banner "Step 1: fitVAR.R chain $SLURM_ARRAY_TASK_ID"
  out="$STEP1_DIR/dsem_step1_${SLURM_ARRAY_TASK_ID}.rds"
  if [[ -s "$out" ]]; then echo "Exists, skipping: $out"; exit 0; fi
  load_r; use_private_tmp
  cd "$STEP1_DIR"
  # Task id is used as both the Stan seed and the output-file suffix.
  Rscript "$REPO_DIR/fitVAR.R" "$SLURM_ARRAY_TASK_ID"
  ;;

granger)
  banner "Step 1b: granger_analysis.R shard $SLURM_ARRAY_TASK_ID/$GRANGER_SHARDS"
  load_r; use_private_tmp
  # granger_analysis.R reads every *.rds in its input folder and writes its
  # own .rds output there, so each shard gets a private folder of symlinks
  # to the Step 1 files only.
  shard_dir="$GRANGER_DIR/shard_$(printf '%03d' "$SLURM_ARRAY_TASK_ID")"
  mkdir -p "$shard_dir"
  find "$shard_dir" -maxdepth 1 -type l -name 'dsem_step1_*.rds' -delete
  ln -s "$STEP1_DIR"/dsem_step1_*.rds "$shard_dir"/
  echo "Using $(find "$shard_dir" -maxdepth 1 -name 'dsem_step1_*.rds' | wc -l) Step 1 files"
  Rscript "$REPO_DIR/granger_analysis.R" "$shard_dir" "$GRANGER_PREFIX" \
    array_task_id="$SLURM_ARRAY_TASK_ID" array_task_count="$GRANGER_SHARDS" \
    block_rows="$GRANGER_BLOCK_ROWS"
  ;;

granger_combine)
  banner "Step 1b: combine Granger shards"
  load_r
  Rscript - "$GRANGER_DIR" "$GRANGER_PREFIX" <<'EOF'
args   <- commandArgs(trailingOnly = TRUE)
gdir   <- args[1]; prefix <- args[2]
files  <- list.files(gdir, pattern = paste0("^", prefix, ".*_pairwise_elpd\\.csv$"),
                     recursive = TRUE, full.names = TRUE)
files  <- files[dirname(files) != gdir]          # shard outputs only
if (length(files) == 0) stop("No Granger shard outputs found under ", gdir)
res <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
res <- res[order(res$elpd_diff_full_minus_reduced, decreasing = TRUE), ]
out <- file.path(gdir, paste0(prefix, "_pairwise_elpd.csv"))
write.csv(res, out, row.names = FALSE)
saveRDS(res, sub("\\.csv$", ".rds", out))
cat("Combined", length(files), "shard files ->", nrow(res), "directed pairs (paper: 72)\n")
pos  <- res$supports_x_granger_causing_y & res$mean_coef_x_to_y > 0
cat("ELPD improved AND positive mean coefficient:", sum(pos), "(paper: 19)\n")
cat("  ...of which z > 2:", sum(pos & res$strong_support_diff_gt_2se), "(paper: 10)\n")
for (z in intersect(c("z_child", "z_block"), names(res)))
  cat("  ...of which", z, "> 2:", sum(pos & res[[z]] > 2, na.rm = TRUE), "\n")
cat("Strongest edge:", res$predictor_x[1], "->", res$outcome_y[1], "(paper: gaze -> contingent)\n")
cat("Saved", out, "\n")
EOF
  ;;

step2a)
  b="${BEHAVIOURS[$SLURM_ARRAY_TASK_ID]}"
  [[ -n "$b" ]] || { echo "No behaviour for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
  banner "Step 2a: select_terms_for_individual_behaviours.R ($b)"
  load_r; use_private_tmp
  cd "$REPO_DIR"                 # script reads ./CDI_data.csv
  Rscript "$REPO_DIR/select_terms_for_individual_behaviours.R" "$b" "$STEP1_DIR" \
    ${STEP2_SEED:+seed=$STEP2_SEED}
  ;;

convergence)
  banner "Convergence check of Step 1 fits"
  load_r; use_private_tmp
  Rscript "$REPO_DIR/check_convergence.R" "$STEP1_DIR" "$RUN_DIR/CONVERGENCE"
  ;;

granger_direction)
  banner "Step 1c: granger_direction.R"
  load_r
  Rscript "$REPO_DIR/granger_direction.R" "$RUN_DIR"
  ;;

figures)
  banner "Figures 1 and 2"
  load_r
  Rscript "$REPO_DIR/make_figures.R" "$RUN_DIR" "$RUN_DIR/FIGURES" figure="$FIGURE" \
    rule="$FIG1_RULE" prob="$FIG1_PROB" z="$FIG1_Z" show="$FIG1_SHOW" threshold="$FIG2_THRESHOLD" \
    weight="$FIG2_WEIGHT" min_gain="$FIG2_MIN_GAIN"
  ;;

tables)
  banner "Supplementary tables"
  load_r
  Rscript "$REPO_DIR/make_tables.R" "$RUN_DIR" "$RUN_DIR/TABLES"
  ;;

bridge)
  banner "make_bridge.R (Step 1 summaries for the permutation runs)"
  load_r
  mkdir -p "$PERM_DIR"
  Rscript "$REPO_DIR/make_bridge.R" "$STEP1_DIR" "$BRIDGE_FILE"
  ;;

perm2a)
  n_beh=${#BEHAVIOURS[@]}
  k=$(( SLURM_ARRAY_TASK_ID / n_beh + ${PERM_START:-1} ))
  b="${BEHAVIOURS[$(( SLURM_ARRAY_TASK_ID % n_beh ))]}"
  pdir="$PERM_DIR/perm_$(printf '%03d' "$k")"
  banner "Permutation $k, Step 2a ($b)"
  out="$pdir/dsem_selection_${b}_single_behaviour.rds"
  if [[ -s "$out" ]]; then echo "Exists, skipping: $out"; exit 0; fi
  [[ -s "$BRIDGE_FILE" ]] || { echo "Missing $BRIDGE_FILE (run the bridge stage)" >&2; exit 1; }
  load_r; use_private_tmp
  mkdir -p "$pdir"
  cd "$REPO_DIR"                 # script reads ./CDI_data.csv
  Rscript "$REPO_DIR/select_terms_for_individual_behaviours.R" "$b" "$STEP1_DIR" \
    bridge="$BRIDGE_FILE" out_dir="$pdir" permute_seed=$(( PERM_SEED_BASE + k )) \
    seed=$(( PERM_SEED_BASE + 1000 * k )) chains="$PERM_CHAINS" iter="$PERM_ITER"
  ;;

perm2b)
  k="$SLURM_ARRAY_TASK_ID"
  pdir="$PERM_DIR/perm_$(printf '%03d' "$k")"
  banner "Permutation $k, Step 2b"
  if [[ -s "$pdir/$STEP2B_PREFIX.rds" ]]; then echo "Exists, skipping: $pdir/$STEP2B_PREFIX.rds"; exit 0; fi
  n_sel=$(find "$pdir" -maxdepth 1 -name 'dsem_selection_*_single_behaviour.rds' 2>/dev/null | wc -l)
  if [[ "$n_sel" -ne ${#BEHAVIOURS[@]} ]]; then
    echo "Only $n_sel of ${#BEHAVIOURS[@]} Step 2a files in $pdir; skipping (resubmit perm2a, then perm2b)" >&2; exit 1
  fi
  load_r; use_private_tmp
  cd "$REPO_DIR"
  Rscript "$REPO_DIR/final_term_weighting.R" "$pdir" "$STEP1_DIR" "$STEP2B_PREFIX" \
    bridge="$BRIDGE_FILE" permute_seed=$(( PERM_SEED_BASE + k )) seed=$(( PERM_SEED_BASE + 1000 * k )) \
    blocks="$STEP2B_BLOCKS" min_inner="$STEP2B_MIN_INNER" weights="$STEP2_WEIGHTS" \
    chains="$PERM_CHAINS" iter="$PERM_ITER"
  ;;

sb_fit)
  banner "single_block_test.R run $SLURM_ARRAY_TASK_ID (0 = observed)"
  load_r; use_private_tmp
  cd "$REPO_DIR"                 # script reads ./CDI_data.csv
  Rscript "$REPO_DIR/single_block_test.R" fit "$RUN_DIR" "$SLURM_ARRAY_TASK_ID" \
    perm_seed_base="$PERM_SEED_BASE" chains=4 iter=500
  ;;

sb_summary)
  banner "single_block_test.R summary"
  load_r
  Rscript "$REPO_DIR/single_block_test.R" summary "$RUN_DIR"
  ;;

perm2b_only)
  k="$SLURM_ARRAY_TASK_ID"
  if [[ "$k" == 0 ]]; then
    pdir="$STEP1_DIR"; perm_args=""; chains_k=4      # observed data, as many chains as Step 2b
  else
    pdir="$PERM_DIR/perm_$(printf '%03d' "$k")"; perm_args="permute_seed=$(( PERM_SEED_BASE + k ))"; chains_k="$PERM_CHAINS"
  fi
  banner "Run $k (0 = observed), Step 2b only (group level)${VARIANT_TAG:+ [$VARIANT_TAG]}"
  if [[ -s "$pdir/$VARIANT_PREFIX.rds" ]]; then echo "Exists, skipping: $pdir/$VARIANT_PREFIX.rds"; exit 0; fi
  load_r; use_private_tmp
  mkdir -p "$pdir"
  cd "$REPO_DIR"
  Rscript "$REPO_DIR/final_term_weighting.R" "$pdir" "$STEP1_DIR" "$VARIANT_PREFIX" \
    bridge="$BRIDGE_FILE" $perm_args seed=$(( PERM_SEED_BASE + 1000 * k )) \
    blocks=all no_selection=1 weights="$STEP2_WEIGHTS" chains="$chains_k" iter="$PERM_ITER" \
    grouping="$STEP2B_GROUPING" ${STEP2B_EXCLUDE:+exclude=$STEP2B_EXCLUDE}
  ;;

perm_summary)
  banner "permutation_null.R"
  load_r
  Rscript "$REPO_DIR/permutation_null.R" "$RUN_DIR" "$VARIANT_PREFIX"
  ;;

step2b)
  banner "Step 2b: final_term_weighting.R"
  load_r; use_private_tmp
  n_sel=$(find "$STEP1_DIR" -maxdepth 1 -name 'dsem_selection_*_single_behaviour.rds' | wc -l)
  echo "Found $n_sel Step 2a selection files (expect ${#BEHAVIOURS[@]})"
  cd "$REPO_DIR"                 # script reads ./CDI_data.csv
  Rscript "$REPO_DIR/final_term_weighting.R" "$STEP1_DIR" "$STEP1_DIR" "$STEP2B_PREFIX" \
    ${STEP2_SEED:+seed=$STEP2_SEED} blocks="$STEP2B_BLOCKS" min_inner="$STEP2B_MIN_INNER" weights="$STEP2_WEIGHTS"
  ;;

*)
  echo "Unknown stage '$STAGE'. Use: submit | step1 | convergence | granger | granger_combine | granger_direction | step2a | step2b | figures | tables | permutations | bridge | perm2a | perm2b | perm_summary | single_block | sb_summary | perm_groups | group_tests" >&2
  exit 1
  ;;
esac

banner "Finished $STAGE"
