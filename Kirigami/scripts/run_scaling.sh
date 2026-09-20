#!/bin/zsh
# WP7b -- drives Kirigami/apps/exp_scaling.jl one (cell, routine) per process, strictly
# sequentially, so that the peak resident set at exit (Sys.maxrss) is that routine's own.
#
# Each invocation is wrapped in the repo's `perl -e alarm` cap (no GNU coreutils `timeout`
# on this machine). A routine that is killed writes no row; the script records the miss in
# capped.csv, and STOPS running that routine at larger sizes for that family -- the cost
# grows monotonically, so continuing only burns the cap again.
#
#   Kirigami/scripts/run_scaling.sh [outdir]      (run from the repo root)
set -u
export PATH=$HOME/.juliaup/bin:$PATH
OUT=${1:-results/scaling}
BIN=(julia --project=Kirigami Kirigami/apps/exp_scaling.jl)
CSV=$OUT/scaling.csv
CAP_CSV=$OUT/capped.csv
mkdir -p $OUT
[[ -f $CAP_CSV ]] || echo "kind,sites,seed,routine,cap_s,outcome" > $CAP_CSV

TARGETS=(50 100 200 500 1000 2000 5000)
SEEDS=(0 1 2)
# routine -> wall cap, seconds. `solve` carries the spec's ~600 s dense-SVD budget.
typeset -A CAP
CAP[solve]=600
CAP[sparse]=1800
CAP[characterize]=1800
CAP[rangemax]=1800

run_stage() {
  local kind=$1 sites=$2 seed=$3 stage=$4 cap=$5
  local t0=$SECONDS
  /usr/bin/perl -e 'alarm shift; exec @ARGV' $cap \
      $BIN[@] --kind $kind --sites $sites --seed $seed --stage $stage --out $OUT --csv $CSV \
      >> $OUT/run.log 2>&1
  local rc=$?
  local dt=$((SECONDS - t0))
  if [[ $rc -ne 0 ]]; then
    echo "$kind,$sites,$seed,$stage,$cap,capped_or_failed_rc${rc}_after_${dt}s" >> $CAP_CSV
    echo "[scaling] CAP  $kind sites=$sites seed=$seed $stage rc=$rc after ${dt}s"
    return 1
  fi
  echo "[scaling] ok   $kind sites=$sites seed=$seed $stage in ${dt}s"
  return 0
}

for kind in delaunay voronoi; do
  # once a routine caps at every seed of one size, it is not tried at larger sizes.
  typeset -A dead
  for r in solve sparse characterize rangemax; do dead[$r]=0; done
  for t in $TARGETS; do
    if [[ $kind == delaunay ]]; then sites=$((t / 2)); else sites=$t; fi
    for stage in solve sparse characterize rangemax; do
      if [[ ${dead[$stage]} -eq 1 ]]; then
        for s in $SEEDS; do
          echo "$kind,$sites,$s,$stage,${CAP[$stage]},skipped_after_cap_at_smaller_size" >> $CAP_CSV
        done
        continue
      fi
      fails=0
      for s in $SEEDS; do
        run_stage $kind $sites $s $stage ${CAP[$stage]} || fails=$((fails + 1))
      done
      if [[ $fails -eq ${#SEEDS[@]} ]]; then
        dead[$stage]=1
        echo "[scaling] routine $stage retired for $kind at sites=$sites (all seeds capped)"
        # characterize and rangemax both need the dense null space; if solve is gone so
        # are they.
        if [[ $stage == solve ]]; then dead[characterize]=1; dead[rangemax]=1; fi
      fi
    done
  done
done
echo "[scaling] done -> $CSV"
