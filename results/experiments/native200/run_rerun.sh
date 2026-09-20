#!/bin/bash
# Native200 WP5 -- the 3600 s rerun of the 362 timed-out + 27 never-dispatched cells.
#
# 389 cells, sorted by |F| ascending and dealt round-robin into 6 shards, so every
# worker reports its cheap cells first.  12-way parallel from 21:50 (orchestrator widened it; the machine is shared; the
# crash-fixed rerun adds 2 more processes, total <= 8).  Each shard appends to its
# own CSV and is restartable: rerunning this script resumes from those CSVs and
# skips every cell already `completed` or `crashed`.
#
#   nohup bash results/experiments/native200/run_rerun.sh > results/experiments/native200/rerun3600/run_all.log 2>&1 &
#
set -u
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
OUT=results/experiments/native200/rerun3600
BIN="julia --project=Kirigami Kirigami/apps/exp_native200.jl"
NSHARD=12
TAG=rerun3600
TIMEOUT=3600

mkdir -p "$OUT/logs" "$OUT/cells"

# ---- progress writer: recomputed from the shard CSVs, so it tracks every cell ---
progress() {
  local total done comp to cr last
  total=$(wc -l < "$OUT/cells_rerun_all.txt" | tr -d ' ')
  {
    echo "native200 rerun3600 -- progress at $(date '+%Y-%m-%d %H:%M:%S')"
    echo "cells in this run: $total (362 timed_out + 27 never dispatched), timeout ${TIMEOUT}s, ${NSHARD}-way"
    if ls "$OUT"/shard_*.csv >/dev/null 2>&1; then
      cat "$OUT"/shard_*.csv "$OUT"/done_6way/*.csv | grep -v '^id,' > /tmp/native200_progress_$$.csv
      done=$(wc -l < /tmp/native200_progress_$$.csv | tr -d ' ')
      comp=$(awk -F, '$7=="completed"' /tmp/native200_progress_$$.csv | wc -l | tr -d ' ')
      to=$(awk -F, '$7=="timed_out"' /tmp/native200_progress_$$.csv | wc -l | tr -d ' ')
      cr=$(awk -F, '$7=="crashed"' /tmp/native200_progress_$$.csv | wc -l | tr -d ' ')
      echo "done $done / $total   completed $comp   timed_out $to   crashed $cr"
      for s in $(seq 0 $((NSHARD-1))); do
        if [ -f "$OUT/shard_$s.csv" ]; then
          local n w
          n=$(grep -vc '^id,' "$OUT/shard_$s.csv")
          w=$(wc -l < "$OUT/cells/shard_$s.txt" | tr -d ' ')
          last=$(grep -v '^id,' "$OUT/shard_$s.csv" | tail -1 | awk -F, '{print "id="$1" "$3" F="$5" -> "$7" ("$17"s)"}')
          echo "  shard $s: $n/$w   last: $last"
        else
          echo "  shard $s: 0/$(wc -l < "$OUT/cells/shard_$s.txt" | tr -d ' ')   last: (none yet)"
        fi
      done
      last=$(cat "$OUT"/shard_*.csv "$OUT"/done_6way/*.csv | grep -v '^id,' | tail -1 | awk -F, '{print "id="$1" kind="$2" variant="$3" F="$5" status="$7" secs="$17}')
      echo "last cell finished (any shard): $last"
      rm -f /tmp/native200_progress_$$.csv
    else
      echo "done 0 / $total   (no shard CSV yet)"
    fi
  } > "$OUT/progress.txt.tmp"
  mv "$OUT/progress.txt.tmp" "$OUT/progress.txt"
}

progress
( while true; do sleep 10; progress; done ) &
PROGRESS_PID=$!
trap 'kill $PROGRESS_PID 2>/dev/null' EXIT

for s in $(seq 0 $((NSHARD-1))); do
  $BIN \
    --cells "$OUT/cells/shard_$s.txt" \
    --resume "$OUT/shard_$s.csv" \
    --resume "$OUT/done_6way/shard_0.csv" --resume "$OUT/done_6way/shard_1.csv" --resume "$OUT/done_6way/shard_2.csv" --resume "$OUT/done_6way/shard_3.csv" --resume "$OUT/done_6way/shard_4.csv" --resume "$OUT/done_6way/shard_5.csv" \
    --csv "$OUT/shard_$s.csv" \
    --tag "$TAG" \
    --timeout "$TIMEOUT" \
    --n 200 --maxf 800 \
    --sigma results/experiments/k5/sigma \
    --cli baseline/native/build/tuttekiri_cli \
    --work "/tmp/kiri_native200_rerun/shard_$s" \
    --out "$OUT" \
    --logdir "$OUT/logs" \
    > "$OUT/logs/run_shard_$s.log" 2>&1 &
done
wait
kill $PROGRESS_PID 2>/dev/null
progress
echo "[run_rerun] all shards finished at $(date)"
