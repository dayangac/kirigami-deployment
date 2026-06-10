#!/bin/bash
# Native200 WP5 -- the 3600 s rerun of the 39 cells that segfaulted the authors'
# CLI, against the F24-patched binary baseline/native/build_fixed/tuttekiri_cli
# (baseline/native/src_fixed/*.cpp = upstream + crashfix.patch; the original
# build/tuttekiri_cli is untouched).
#
# 2-way parallel, so that with run_rerun.sh's 6 workers the total stays at 8.
#
#   nohup bash results/kill/native200/run_crashfix.sh > results/kill/native200/crashfix3600/run_all.log 2>&1 &
#
set -u
ROOT=/Users/emredayangac/Documents/kirigami-experiments
cd "$ROOT"
OUT=results/kill/native200/crashfix3600
BIN=code/build/kill_native200
CLI=baseline/native/build_fixed/tuttekiri_cli
NSHARD=2
TAG=crashfix3600
TIMEOUT=3600

mkdir -p "$OUT/logs" "$OUT/cells"

progress() {
  local total done comp to cr last n w
  total=$(wc -l < "$OUT/cells_crashed_all.txt" | tr -d ' ')
  {
    echo "native200 crashfix3600 -- progress at $(date '+%Y-%m-%d %H:%M:%S')"
    echo "cells in this run: $total (the 39 that segfaulted at 600 s), timeout ${TIMEOUT}s, ${NSHARD}-way, CLI $CLI"
    if ls "$OUT"/shard_*.csv >/dev/null 2>&1; then
      cat "$OUT"/shard_*.csv | grep -v '^id,' > /tmp/native200_cfprogress_$$.csv
      done=$(wc -l < /tmp/native200_cfprogress_$$.csv | tr -d ' ')
      comp=$(awk -F, '$7=="completed"' /tmp/native200_cfprogress_$$.csv | wc -l | tr -d ' ')
      to=$(awk -F, '$7=="timed_out"' /tmp/native200_cfprogress_$$.csv | wc -l | tr -d ' ')
      cr=$(awk -F, '$7=="crashed"' /tmp/native200_cfprogress_$$.csv | wc -l | tr -d ' ')
      echo "done $done / $total   completed $comp   timed_out $to   crashed $cr"
      for s in $(seq 0 $((NSHARD-1))); do
        if [ -f "$OUT/shard_$s.csv" ]; then
          n=$(grep -vc '^id,' "$OUT/shard_$s.csv")
          w=$(wc -l < "$OUT/cells/shard_$s.txt" | tr -d ' ')
          last=$(grep -v '^id,' "$OUT/shard_$s.csv" | tail -1 | awk -F, '{print "id="$1" "$3" F="$5" -> "$7" ("$17"s)"}')
          echo "  shard $s: $n/$w   last: $last"
        else
          echo "  shard $s: 0/$(wc -l < "$OUT/cells/shard_$s.txt" | tr -d ' ')   last: (none yet)"
        fi
      done
      last=$(cat "$OUT"/shard_*.csv | grep -v '^id,' | tail -1 | awk -F, '{print "id="$1" kind="$2" variant="$3" F="$5" status="$7" secs="$17}')
      echo "last cell finished (any shard): $last"
      rm -f /tmp/native200_cfprogress_$$.csv
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
  "$BIN" \
    --cells "$OUT/cells/shard_$s.txt" \
    --resume "$OUT/shard_$s.csv" \
    --csv "$OUT/shard_$s.csv" \
    --tag "$TAG" \
    --timeout "$TIMEOUT" \
    --n 200 --maxf 800 \
    --sigma results/kill/k5/sigma \
    --cli "$CLI" \
    --work "/tmp/kiri_native200_crashfix/shard_$s" \
    --out "$OUT" \
    --logdir "$OUT/logs" \
    > "$OUT/logs/run_shard_$s.log" 2>&1 &
done
wait
kill $PROGRESS_PID 2>/dev/null
progress
echo "[run_crashfix] all shards finished at $(date)"
