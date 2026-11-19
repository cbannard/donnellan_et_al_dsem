#!/bin/bash
# Joint test of SETS of Step 2 predictor blocks against the intercept-only
# model, with the same vocabulary shuffles as the single-block test.
#
#   RUN_DIR=/scratch/.../dsem_session sbatch --array=0-10 run_block_sets.slurm
#
# Task 0 fits the observed vocabulary; task i (1..10) fits shuffles
# (i-1)*100+1 .. i*100. Afterwards, summarise with
#   Rscript single_block_test.R summary $RUN_DIR sets="$SETS" set_tag=$SET_TAG
#
#SBATCH -p himem
#SBATCH -n 1
#SBATCH -t 0-12
#SBATCH --mem=8G
#SBATCH -J dsem_sets
set -euo pipefail
: "${RUN_DIR:?set RUN_DIR}"
SET_TAG="${SET_TAG:-responses}"
SETS="${SETS:-vocal_responses=CV.contingent+non_CV.contingent;all_responses=CV.contingent+gaze.contingent+give.contingent+index_point.contingent+non_CV.contingent+open_point.contingent+other_gesture.contingent+show.contingent}"
PER_TASK="${PER_TASK:-100}"
module purge; module load "${R_MODULE:-apps/gcc/R/4.4.2}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR="$RUN_DIR/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT
cd "${SLURM_SUBMIT_DIR}"          # the repo: script reads ./CDI_data.csv
k=$SLURM_ARRAY_TASK_ID
if (( k == 0 )); then range=0; else range="$(( (k - 1) * PER_TASK + 1 ))-$(( k * PER_TASK ))"; fi
echo "=== block sets [$SET_TAG] runs $range | $(hostname) | $(date) ==="
Rscript single_block_test.R fit "$RUN_DIR" "$range" sets="$SETS" set_tag="$SET_TAG" chains=4 iter=500
