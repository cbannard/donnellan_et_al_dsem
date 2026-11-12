#!/bin/bash
# =============================================================================
# monitor_progress.sh — how far has each stage of a pipeline run got?
#
# Usage (on a login node):
#   bash monitor_progress.sh [RUN_DIR]          default: $RUN_DIR
#   watch -n 300 bash monitor_progress.sh /scratch/q55362cb/dsem_results_v3
#
# Reads only the run folder and the queue; changes nothing. For each stage it
# shows what has finished and, for tasks still running, how far they are:
#   Step 1     sampler iteration of each chain (from its log)
#   Granger    directed pairs scored (8 per shard)
#   Step 2a    models fitted out of 1,024 (counted from CmdStan output files
#              in the task's private temp folder), with an ETA
#   Step 2b    the script's own progress line (models out of 511, ETA)
#   Perms      permutation files finished out of PERM_N x 9 (2a) and PERM_N
#              (2b), plus the progress of running 2a tasks
# =============================================================================

RUN_DIR="${1:-${RUN_DIR:-}}"
[[ -n "$RUN_DIR" && -d "$RUN_DIR" ]] || { echo "Usage: bash monitor_progress.sh <RUN_DIR>" >&2; exit 1; }
LOG="$RUN_DIR/logs"; TMP="$RUN_DIR/tmp"
BEHAVIOURS=(contingent gaze CV non_CV index_point open_point give show other_gesture)
PERM_N="${PERM_N:-100}"
N_COMBOS_2A=1024   # 1,023 block combinations + intercept-only

# Seconds from a Slurm elapsed time ([D-]HH:MM:SS or MM:SS)
secs() {
  local t="$1" d=0
  [[ "$t" == *-* ]] && { d=${t%%-*}; t=${t#*-}; }
  IFS=: read -r a b c <<< "$t"
  if [[ -z "$c" ]]; then echo $(( d*86400 + 10#$a*60 + 10#$b )); else echo $(( d*86400 + 10#$a*3600 + 10#$b*60 + 10#$c )); fi
}
hm() { printf "%dh%02dm" $(( $1 / 3600 )) $(( ($1 % 3600) / 60 )); }

# Running tasks of this user: jobid_arrayidx name elapsed jobid(raw)
running=$(squeue -u "$USER" -h -t R -o "%i %j %M %A" 2>/dev/null)
pending=$(squeue -u "$USER" -h -t PD -r -o "%j" 2>/dev/null | sort | uniq -c)

echo "Run folder: $RUN_DIR    ($(date '+%a %d %b %H:%M'))"
echo
echo "QUEUE (your jobs)"
squeue -u "$USER" -h -r -o "%j %T" 2>/dev/null | sort | uniq -c | awk '{printf "  %-22s %-9s %5d\n", $2, $3, $1}'
echo

# ---- Step 1 ------------------------------------------------------------------
n1=$(ls "$RUN_DIR"/MODELS-STEP1/dsem_step1_*.rds 2>/dev/null | wc -l)
echo "STEP 1  fits finished: $n1"
for f in "$LOG"/step1_*.out; do
  [[ -f "$f" ]] || continue
  last=$(grep -a "Iteration:" "$f" 2>/dev/null | tail -1)
  done_flag=$(grep -c "Finished step1" "$f" 2>/dev/null)
  [[ -n "$last" && "$done_flag" == 0 ]] && echo "  $(basename "$f" .out): ${last#*Iteration:}"
done
echo

# ---- Granger -----------------------------------------------------------------
if [[ -f "$RUN_DIR/GRANGER/granger_loocv_pairwise_elpd.csv" ]]; then
  echo "GRANGER combined: $(($(wc -l < "$RUN_DIR/GRANGER/granger_loocv_pairwise_elpd.csv") - 1)) pairs"
fi
tot=0; shards=0
for f in "$LOG"/granger_[0-9]*.out; do
  [[ -f "$f" ]] || continue
  n=$(grep -c "^  Pair " "$f" 2>/dev/null); tot=$((tot + n)); shards=$((shards + 1))
  st=$(grep -q "Finished granger" "$f" && echo done || echo "")
  printf "  %-12s %d/8 pairs %s\n" "$(basename "$f" .out)" "$n" "$st"
done
[[ $shards -gt 0 ]] && echo "  GRANGER total: $tot / 72 pairs scored"
echo

# ---- Step 2a / perm2a: count fitted models in a task's temp folder ------------
# cmdstanr writes one CSV per chain per model into the task's temp folder.
fits_done() {   # $1 = raw job id, $2 = array index, $3 = chains per model
  local d="$TMP/${1}_${2}"
  [[ -d "$d" ]] || { echo 0; return; }
  local n; n=$(find "$d" -name '*.csv' 2>/dev/null | wc -l)
  echo $(( n / $3 ))
}
show_2a_task() {   # $1 label, $2 raw id, $3 idx, $4 elapsed, $5 chains
  local f; f=$(fits_done "$2" "$3" "$5")
  local s; s=$(secs "$4")
  local eta="?"
  if [[ $f -gt 5 ]]; then eta=$(hm $(( s * (N_COMBOS_2A - f) / f ))); fi
  printf "  %-26s %4d / %d models (%3d%%)  elapsed %-8s ETA %s\n" "$1" "$f" "$N_COMBOS_2A" \
    $(( 100 * f / N_COMBOS_2A )) "$(hm "$s")" "$eta"
}

n2a=$(ls "$RUN_DIR"/MODELS-STEP1/dsem_selection_*_single_behaviour.rds 2>/dev/null | wc -l)
echo "STEP 2a  behaviours finished: $n2a / 9"
while read -r jid name el raw; do
  [[ "$name" == dsem_step2a ]] || continue
  idx=${jid##*_}
  show_2a_task "${BEHAVIOURS[$idx]}" "$raw" "$idx" "$el" 4
done <<< "$running"
echo

# ---- Step 2b -----------------------------------------------------------------
echo "STEP 2b"
if [[ -f "$LOG/step2b.out" ]]; then
  if grep -q "Finished step2b" "$LOG/step2b.out"; then echo "  finished"
  else p=$(grep -a "^Progress" "$LOG/step2b.out" | tail -1); echo "  ${p:-started, no progress line yet}"; fi
else echo "  not started"; fi
echo

# ---- Permutations ------------------------------------------------------------
PD="$RUN_DIR/PERMUTATIONS"
if [[ -d "$PD" ]]; then
  echo "PERMUTATIONS (PERM_N=$PERM_N; set PERM_N=... if different)"
  [[ -s "$PD/bridge.rds" ]] && echo "  bridge: done" || echo "  bridge: not done"
  p2a=$(ls "$PD"/perm_*/dsem_selection_*_single_behaviour.rds 2>/dev/null | wc -l)
  p2b=$(ls "$PD"/perm_*/dsem_selection_from_single_behaviour_sets.rds 2>/dev/null | wc -l)
  printf "  perm 2a tasks finished: %d / %d (%d%%)\n" "$p2a" $((PERM_N * 9)) $(( 100 * p2a / (PERM_N * 9) ))
  printf "  perm 2b runs finished:  %d / %d\n" "$p2b" "$PERM_N"
  # Average run time of finished perm2a tasks (from Slurm accounting;
  # tasks under 10 minutes were skipped because their output already existed)
  mapfile -t durs < <(sacct -u "$USER" -X -n --name=dsem_perm2a -s CD -o ElapsedRaw \
                        -S "$(date -d '-30 days' +%Y-%m-%d)" 2>/dev/null | awk '$1 > 600 {print $1}')
  nrun=$(awk '$2=="dsem_perm2a"' <<< "$running" | wc -l)
  if [[ ${#durs[@]} -gt 0 ]]; then
    avg=$(( $(IFS=+; echo "$(( ${durs[*]} ))") / ${#durs[@]} ))
    left=$(( PERM_N * 9 - p2a ))
    par=$(( nrun > 0 ? nrun : 1 ))
    echo "  average perm 2a task: $(hm $avg) over ${#durs[@]} tasks;" \
         "rough time to finish 2a at $par running: $(hm $(( avg * left / par )))"
  fi
  [[ -f "$PD/permutation_null_report.txt" ]] && echo "  summary: done ($PD/permutation_null_report.txt)"
  # Running perm2a tasks
  if [[ $nrun -gt 0 ]]; then
    echo "  running perm 2a tasks ($nrun):"
    while read -r jid name el raw; do
      [[ "$name" == dsem_perm2a ]] || continue
      idx=${jid##*_}; k=$(( idx / 9 + 1 ))
      show_2a_task "perm $k ${BEHAVIOURS[$(( idx % 9 ))]}" "$raw" "$idx" "$el" "${PERM_CHAINS:-2}"
    done <<< "$running" | head -20
    [[ $nrun -gt 20 ]] && echo "  ... ($((nrun - 20)) more)"
  fi
fi

# ---- Failures ----------------------------------------------------------------
fails=$(sacct -u "$USER" -X -n -S "$(date -d '-2 days' +%Y-%m-%d)" -s F,OOM,TO,NF,CA \
          -o JobID%20,JobName%20,State%14 2>/dev/null | grep dsem_)
if [[ -n "$fails" ]]; then
  echo
  echo "FAILED / CANCELLED in the last 2 days (check \$RUN_DIR/logs before resubmitting):"
  echo "$fails" | sort -k2 | awk '{print "  "$0}' | tail -30
fi
