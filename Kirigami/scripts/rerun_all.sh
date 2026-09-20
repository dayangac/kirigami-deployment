#!/bin/bash
# rerun_all.sh -- rerun every experiment driver on its FULL population, with the
# arguments the committed results were produced with (results/experiments/EXPERIMENTS.md
# "Reproducing" blocks, ../paper/techreport/techreport.tex (outside the repository) Sec. "Reproduction",
# Kirigami/scripts/run_e1.sh / run_scaling.sh, results/*/run.log), cheapest first.
# Outputs land in the apps' default results/ directories (results/experiments/<x>/,
# results/final/e1/, results/{regime,scaling,yield}/), overwriting the committed copies,
# so run it on a clean tree and compare against the committed versions with git.
# Optimiser endpoints (K6, B4, K9c, K7 C3, K1c's sweep) are path dependent across
# floating-point environments (docs/NUMERICS.md); the combinatorial columns and the
# referee are not.
#
#   nohup bash Kirigami/scripts/rerun_all.sh > results/logs/rerun_all.log 2>&1 &
#   JULIA=/path/to/arm64/julia NSHARDS=12 bash Kirigami/scripts/rerun_all.sh [--only k5,k6] [--dry-run]
#
# Requires an arm64 Julia (bit-exact trig against the frozen corpora, docs/NUMERICS.md;
# an x86_64 Julia under Rosetta runs everything but tilings/relaxation may round 1 ulp
# differently, see Kirigami/src/core/generators.jl). Estimated wall times are the per-row `secs` sums
# of the committed CSVs (results/*), divided by the shard count where the run is sharded;
# native200 and K2b call the authors' native CLI (baseline/native/build/tuttekiri_cli) and
# run last.
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
JULIA=${JULIA:-julia}
NSHARDS=${NSHARDS:-12}
LOGS=results/logs
ONLY=""
DRY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --only=*) ONLY="${1#--only=}" ;;
    --only) shift; ONLY="$1" ;;
    --dry-run) DRY=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
mkdir -p "$LOGS"

# ---- toolchain check ---------------------------------------------------------------
ARCH=$("$JULIA" -e 'print(Sys.ARCH)' 2>/dev/null || echo unknown)
if [ "$ARCH" != "aarch64" ]; then
  echo "WARNING: $JULIA is $ARCH, not arm64 -- the corpora are arm64 libm; tiling/relaxation trig may differ by 1 ulp." >&2
  echo "         Set JULIA=<juliaup arm64 julia> for the bit-exact rerun." >&2
fi
APP="$JULIA --project=Kirigami Kirigami/apps"
SCR="$JULIA --project=Kirigami/scripts Kirigami/scripts"

# ---- helpers --------------------------------------------------------------------------
selected() {  # selected NAME -> 0 if NAME (or its prefix before `_`, e.g. k7 for k7_main) is in --only
  [ -z "$ONLY" ] && return 0
  case ",$ONLY," in *",$1,"*) return 0 ;; esac
  case ",$ONLY," in *",${1%%_*},"*) return 0 ;; esac
  return 1
}
stamp() { date '+%Y-%m-%d %H:%M:%S'; }
has_driver() { ls Kirigami/apps/exp_"$1".jl Kirigami/apps/exp_"$1"_*.jl >/dev/null 2>&1; }  # exp_<tag>[_<topic>].jl
run_app() {  # run_app NAME EST "command ..." : one nohup'd process, its own log, wait
  local name=$1 est=$2; shift 2
  selected "$name" || return 0
  if ! has_driver "${name%%[_/]*}" && ! has_driver "$name"; then
    echo "[$(stamp)] SKIP $name: no Julia driver yet ($*)"; return 0
  fi
  echo "[$(stamp)] START $name  (est. wall $est)  :: $*"
  [ "$DRY" = 1 ] && return 0
  local t0=$SECONDS
  nohup bash -c "$*" > "$LOGS/$name.log" 2>&1
  local rc=$?
  echo "[$(stamp)] END   $name  rc=$rc  wall $((SECONDS - t0)) s  -> $LOGS/$name.log"
}
run_sharded() {  # run_sharded NAME EST FLAG "command ..." : NSHARDS concurrent processes
  local name=$1 est=$2 flag=$3; shift 3
  selected "$name" || return 0
  if ! has_driver "${name%%[_/]*}"; then
    echo "[$(stamp)] SKIP $name: no Julia driver yet ($*)"; return 0
  fi
  echo "[$(stamp)] START $name  ($NSHARDS shards, est. wall $est)  :: $* --shard I $flag $NSHARDS"
  [ "$DRY" = 1 ] && return 0
  local t0=$SECONDS
  for i in $(seq 0 $((NSHARDS - 1))); do
    nohup bash -c "$* --shard $i $flag $NSHARDS" > "$LOGS/${name}_shard_$i.log" 2>&1 &
  done
  wait
  echo "[$(stamp)] END   $name  wall $((SECONDS - t0)) s  -> $LOGS/${name}_shard_*.log"
}
merge_csv() {  # merge_csv OUT first.csv others...  (one header, then every data row)
  local out=$1; shift
  [ "$DRY" = 1 ] && return 0
  head -n1 "$1" > "$out"
  for f in "$@"; do tail -n +2 "$f" >> "$out"; done
  echo "        merged $(wc -l < "$out" | tr -d ' ') lines -> $out"
}

echo "rerun_all: $(stamp)  julia=$JULIA ($ARCH)  nshards=$NSHARDS  root=$ROOT"

# ======================================================================================
# Cheapest first. Dependencies: k1a -> {k1b, k2c, k8a} (shape cache results/experiments/cache,
# shape-cache layout); k5 -> k6 (sigma_def, default = corpus) -> {t1, k9, k9b,
# k9c, b4, k8a_recheck} (K6 cache); k9 + k9b -> k9c (their CSVs). Every app defaults to the
# frozen corpus populations, so the order below only matters for the caches.
# ======================================================================================

# --- seconds-scale ---------------------------------------------------------------------
run_app f23      "~1 s"      "$APP/exp_f23_boundary_row.jl"
run_app k2a      "~7 s"      "$APP/exp_k2a_thetamax_bisection.jl --n 292"
run_sharded jitter "~2 s"    --nshards "$APP/exp_a3_jitter.jl"
if selected jitter && [ "$DRY" = 0 ]; then
  merge_csv results/experiments/jitter/jitter.csv results/experiments/jitter/jitter_[0-9]*.csv
  run_app jitter_analyze "~1 s" "$APP/exp_a3_jitter.jl --analyze results/experiments/jitter/jitter.csv"
fi
run_app k3a_recheck "~8 s (committed run: 8 s)"   "$APP/exp_k3a_core_mobility_recheck.jl"
run_app k1b      "~14 s (committed run: 14 s)"    "$APP/exp_k1b_harmonic_identity.jl --n 200 --shapes 20"
run_app k1a      "~36 s (committed run: 36 s)"    "$APP/exp_k1a_projection_validity.jl --n 200 --samples 10000"
run_app k1c      "~31 s, native CLI"   "$APP/exp_k1c_eq9_false_negatives.jl"
run_app k2c      "~47 s (committed run: 47 s)"    "$APP/exp_k2c_active_set_locality.jl --n 500"
run_app t1       "~7 s off the K6 cache" "$APP/exp_t1_balance_defect.jl"
run_app b3       "~1 min"              "$APP/exp_b3_expansion_budget.jl --stage run"

# --- minutes-scale ----------------------------------------------------------------------
run_app k6       "~2 min (committed run: 96 s + native)"  "$APP/exp_k6_zero_plus_repair.jl --n 200 --maxf 800"
run_app k5       "~5 min (committed run: 307 s)"          "$APP/exp_k5_orientation.jl --n 200 --maxf 800"
run_sharded e1   "~6 min / $NSHARDS (committed run: 352 s total)" --nshards "$APP/exp_e1_thetamax_validation.jl --n-random 900"
if selected e1 && [ "$DRY" = 0 ]; then
  merge_csv results/final/e1/e1.csv results/final/e1/e1_[0-9]*.csv
fi
run_sharded k8a  "~11 min / $NSHARDS (committed run: 658 s total)" --nshard "$APP/exp_k8a_expansive_cone.jl --n 100 --x0 --dual-iters 100000"
if selected k8a && [ "$DRY" = 0 ]; then
  merge_csv results/experiments/k8a/k8a.csv results/experiments/k8a/k8a_[0-9]*.csv
  run_app k8a_summary "~1 s" "$SCR/summarise_k8a.jl results/experiments/k8a"
fi
run_app k8a_recheck "~5 min"           "$APP/exp_k8a_expansive_cone_recheck.jl"
run_app k3a      "~24 min (committed run: 1410 s)" "$APP/exp_k3a_core_mobility.jl --n 500"
# K7: main + c3 + perdbg sharded (wall set by voronoi_torus_11_n200, ~30 min alone), then
# the three single-process stages
run_sharded k7_main   "~30 min (committed run wall, 12 shards)" --nshard "$APP/exp_k7_periodic_jacobian.jl --stage main"
run_sharded k7_c3     "~30 min"                      --nshard "$APP/exp_k7_periodic_jacobian.jl --stage c3"
run_sharded k7_perdbg "~5 min"                       --nshard "$APP/exp_k7_periodic_jacobian.jl --stage perdbg"
if selected k7_main && [ "$DRY" = 0 ]; then
  merge_csv results/experiments/k7/k7_main.csv results/experiments/k7/k7_main_[0-9]*.csv
  merge_csv results/experiments/k7/k7_c3_all.csv results/experiments/k7/k7_c3_[0-9]*.csv
  merge_csv results/experiments/k7/k7_periodicity.csv results/experiments/k7/k7_periodicity_[0-9]*.csv
fi
run_app k7_bounded "~1 min" "$APP/exp_k7_periodic_jacobian.jl --stage bounded --out results/experiments/k7"
run_app k7_nu      "~1 min" "$APP/exp_k7_periodic_jacobian.jl --stage nu --out results/experiments/k7"
run_app k7_c4sweep "~1 min" "$APP/exp_k7_periodic_jacobian.jl --stage c4sweep --out results/experiments/k7"

# --- hours-scale ------------------------------------------------------------------------
run_sharded b4   "~72 min / $NSHARDS (committed run: 4313 s total)" --nshards "$APP/exp_b4_boundary_vs_graph.jl --n 200"
if selected b4 && [ "$DRY" = 0 ]; then
  merge_csv results/experiments/b4/b4.csv results/experiments/b4/shard_[0-9]*.csv
fi
run_sharded k9   "~96 min / $NSHARDS (committed run: 5774 s total)" --nshards "$APP/exp_k9_convex_embedding.jl --n 200"
run_app k9_aggregate "~1 s" "$APP/exp_k9_convex_embedding.jl --aggregate --nshards $NSHARDS"
run_sharded k9b  "~7.7 h / $NSHARDS (committed run: 27685 s total)" --nshards "$APP/exp_k9b_convex_embedding_search.jl --n 200"
run_app k9b_aggregate "~1 s" "$APP/exp_k9b_convex_embedding_search.jl --aggregate --nshards $NSHARDS"
run_sharded k9c  "~5.4 h / $NSHARDS (committed run: 19318 s total)" --nshards "$APP/exp_k9c_range_embedding.jl --n 200"
run_app k9c_aggregate "~1 s" "$APP/exp_k9c_range_embedding.jl --aggregate --nshards $NSHARDS"
run_app scaling "hours; sequential cells with wall caps (see run_scaling.sh)" \
  "zsh Kirigami/scripts/run_scaling.sh results/scaling"
# regime: the committed run used 4 shards for the constructive arms and 4 for the native cells
REGIME_SHARDS=4
if selected regime; then
  save=$NSHARDS; NSHARDS=$REGIME_SHARDS
  run_sharded regime_pop    "~1 h / 4 (constructive arms)" --nshards "$APP/exp_regime.jl --mode pop"
  run_app     regime_ref    "~10 min"                       "$APP/exp_regime.jl --mode ref"
  NSHARDS=$save
fi
run_app yield    "hours"  "$APP/exp_yield_features.jl --mode fresh"
run_app basin    "hours"  "$APP/exp_basin.jl"

# --- native CLI last (hours; baseline/native/build/tuttekiri_cli must exist) ------------
run_app k2b      "~26 min (committed run: 1517 s), native CLI" "$APP/exp_k2b_range_margin.jl --maxf 160"
# regime native cells (4 shards, 600 s timeout each), then the merge
if selected regime; then
  save=$NSHARDS; NSHARDS=$REGIME_SHARDS
  run_sharded regime_native "hours, native CLI (600 s timeout per cell)" --nshards "$APP/exp_regime.jl --mode native --timeout 600"
  NSHARDS=$save
  run_app     regime_aggregate "~1 s" "$APP/exp_regime.jl --mode aggregate --merge-shards $REGIME_SHARDS"
fi
run_sharded native200 "~65 h / $NSHARDS CPU (committed run: 234293 s total, 600 s timeout per cell)" --nshards \
  "$APP/exp_native200.jl --n 200 --maxf 800 --timeout 600 --out results/experiments/native200/shards --sigma results/experiments/k5/sigma"

# ---- figures and derived summaries ----------------------------------------------------
if [ "$DRY" = 0 ] && [ -z "$ONLY" ]; then
  echo "[$(stamp)] figures"
  for p in plot_experiments plot_e1 plot_jitter plot_k7 plot_k9b plot_k9c plot_regime plot_scaling plot_final; do
    $SCR/$p.jl > "$LOGS/$p.log" 2>&1 || echo "  $p failed (see $LOGS/$p.log)"
  done
  $SCR/plot_k6.jl results/experiments/k6/k6.csv -o results/experiments/k6/k6_zero_plus.png > "$LOGS/plot_k6.log" 2>&1
  echo "[$(stamp)] done"
fi
