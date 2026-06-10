#!/bin/bash
# TODO(julia-port): invoked C++ binary kill_e1; replace with Kirigami/apps/kill_e1.jl
# Runs E1 (results/final/e1) 12-way sharded, then merges the shard CSVs.
set -euo pipefail
cd "$(dirname "$0")/../build"
OUT=../../results/final/e1
mkdir -p "$OUT/shards"
NSHARDS=12
NRANDOM=${1:-900}

seq 0 $((NSHARDS-1)) | xargs -P "$NSHARDS" -I{} \
  ./kill_e1 --shard {} --nshards $NSHARDS --n-random "$NRANDOM" --out "$OUT/shards"

# merge: one header, then all data rows
head -n1 "$OUT/shards/e1_0.csv" > "$OUT/e1.csv"
for f in "$OUT"/shards/e1_*.csv; do
  tail -n +2 "$f" >> "$OUT/e1.csv"
done
echo "merged: $(wc -l < "$OUT/e1.csv") lines (incl header)"
