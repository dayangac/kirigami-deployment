#!/bin/sh
# Merge the fresh-run shards into results/yield/fresh.csv (header once).
cd "$(dirname "$0")/../.."
head -1 results/yield/fresh_shard_0.csv > results/yield/fresh.csv
for s in 0 1 2 3; do
  [ -f results/yield/fresh_shard_$s.csv ] && awk 'NR>1' results/yield/fresh_shard_$s.csv >> results/yield/fresh.csv
done
echo "fresh.csv rows: $(( $(wc -l < results/yield/fresh.csv) - 1 ))"
