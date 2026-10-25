#!/bin/bash --login
# =============================================================================
# run_dsem_pipeline.slurm
#
# Runs the full Study 1 (DSEM) pipeline from the donnellan_et_al_dsem repo
# (Donnellan et al., Phil. Trans. R. Soc. B, RSTB-2025-0407):
#
#   Step 1   fitVAR.R                                 16 single-chain fits (1000 warmup + 1000 draws each)
#   Step 1b  granger_analysis.R                       Figure 1 / 72 directed pairs
#   Step 2a  select_terms_for_individual_behaviours.R one task per behaviour
#   Step 2b  final_term_weighting.R                   Figure 2 inclusion weights
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
#       step1 (array 1-16)  ──┬── granger (array) ── granger_combine
#                             └── step2a  (array 0-8) ── step2b
#
#   To re-run a single stage (e.g. after fixing a failure):
#
#       bash run_dsem_pipeline.slurm step2a        # also: step1, granger,
#                                               #   granger_combine, step2b
#
#   Step 2 seeds: by default each Step 2 job draws a random seed and records
#   it (log + output .rds). To reproduce a run:  STEP2_SEED=12345 bash run_dsem_pipeline.slurm step2b
#
#   Step 1 skips chains whose dsem_step1_<n>.rds already exists, so
#   re-submitting step1 only reruns missing or failed chains.
#
# OUTPUT (all under $RUN_DIR, default donnellan_et_al_dsem/results/)
#   MODELS-STEP1/dsem_step1_<n>.rds                     Step 1 posterior draws
#   GRANGER/granger_loocv_pairwise_elpd.csv             Step 1b, all 72 pairs
#   MODELS-STEP1/dsem_selection_<b>_single_behaviour.rds  Step 2a
#   MODELS-STEP1/dsem_selection_from_single_behaviour_sets_stacking_coefs.csv
#                                                       Step 2b final weights
#   logs/                                               one log per task
#
# Cluster settings below match the original runs (CSF3-style "himem"
# partition, module apps/gcc/R/4.4.2, 7-day limit). CmdStan (>= 2.33) and the
# R packages cmdstanr, posterior, loo, dplyr, tidyr, readr must already be
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

STEP2_SEED="${STEP2_SEED:-}"         # blank = each Step 2 job draws a random seed
                                      #   (printed in its log and saved in its .rds);
                                      #   set a number to reproduce an earlier run

BEHAVIOURS=(contingent gaze CV non_CV index_point open_point give show other_gesture)
GRANGER_PREFIX="granger_loocv"
STEP2B_PREFIX="dsem_selection_from_single_behaviour_sets"
# -----------------------------------------------------------------------------

set -eo pipefail

STAGE="${STAGE:-${1:-submit}}"

# Paths. When submitting, REPO_DIR defaults to the folder holding this script;
# compute stages receive REPO_DIR/RUN_DIR through --export.
if [[ -z "${REPO_DIR:-}" ]]; then
  REPO_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
fi
RUN_DIR="${RUN_DIR:-$REPO_DIR/results}"
STEP1_DIR="$RUN_DIR/MODELS-STEP1"
GRANGER_DIR="$RUN_DIR/GRANGER"
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
    --export=ALL,STAGE="$stage",REPO_DIR="$REPO_DIR",RUN_DIR="$RUN_DIR",STEP2_SEED="$STEP2_SEED" \
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
           final_term_weighting.R CDI_data.csv; do
    [[ -f "$REPO_DIR/$f" ]] || { echo "Missing $REPO_DIR/$f" >&2; exit 1; }
  done
  if [[ ! -d "$REPO_DIR/data_full" ]]; then
    [[ -f "$REPO_DIR/data_full.zip" ]] || { echo "Missing data_full/ and data_full.zip" >&2; exit 1; }
    echo "Unzipping data_full.zip ..."
    (cd "$REPO_DIR" && unzip -q data_full.zip -x '__MACOSX/*')
  fi
  n_csv=$(find "$REPO_DIR/data_full" -name 'p*_1[12]m.csv' | wc -l)
  echo "Found $n_csv participant CSV files in data_full/ (expect 116 = 58 children x 2 sessions)"

  if command -v module >/dev/null 2>&1; then
    load_r >/dev/null 2>&1 || true
  fi
  if command -v Rscript >/dev/null 2>&1; then
    Rscript -e 'for (p in c("cmdstanr","posterior","loo","dplyr","tidyr","readr")) if (!requireNamespace(p, quietly=TRUE)) stop("R package missing: ", p); v <- cmdstanr::cmdstan_version(); cat("CmdStan", v, "\n"); if (package_version(v) < "2.33") stop("CmdStan >= 2.33 needed for Pathfinder")' \
      || { echo "Preflight R check failed (see above)." >&2; exit 1; }
  else
    echo "WARNING: Rscript not found on this node; skipping package check."
  fi

  mkdir -p "$STEP1_DIR" "$GRANGER_DIR" "$LOG_DIR" "$RUN_DIR/tmp"
  ln -sfn "$REPO_DIR/data_full" "$STEP1_DIR/data_full"   # fitVAR.R reads ./data_full

  n_beh=${#BEHAVIOURS[@]}
  dep() { [[ "$WHAT" == submit && -n "$1" ]] && echo "--dependency=afterok:$1"; }

  J1=""; JG=""; J2A=""
  if [[ "$WHAT" == submit || "$WHAT" == step1 ]]; then
    J1=$(submit_stage step1 -J dsem_step1 -t "$TIME_LONG" $(mem_flag "$MEM_STEP1") \
           --array=1-"$N_CHAINS" -o "$LOG_DIR/step1_%a.out")
    echo "Step 1   (fitVAR, $N_CHAINS chains):           job $J1"
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
    array_task_id="$SLURM_ARRAY_TASK_ID" array_task_count="$GRANGER_SHARDS"
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

step2b)
  banner "Step 2b: final_term_weighting.R"
  load_r; use_private_tmp
  n_sel=$(find "$STEP1_DIR" -maxdepth 1 -name 'dsem_selection_*_single_behaviour.rds' | wc -l)
  echo "Found $n_sel Step 2a selection files (expect ${#BEHAVIOURS[@]})"
  cd "$REPO_DIR"                 # script reads ./CDI_data.csv
  Rscript "$REPO_DIR/final_term_weighting.R" "$STEP1_DIR" "$STEP1_DIR" "$STEP2B_PREFIX" \
    ${STEP2_SEED:+seed=$STEP2_SEED}
  ;;

*)
  echo "Unknown stage '$STAGE'. Use: submit | step1 | granger | granger_combine | step2a | step2b" >&2
  exit 1
  ;;
esac

banner "Finished $STAGE"
