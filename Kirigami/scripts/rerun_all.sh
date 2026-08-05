#!/bin/bash
# rerun_all.sh -- rerun every experiment driver of the Julia port on its FULL population,
# with the arguments the C++ runs used (results/kill/KILL_REPORT.md "Reproducing" blocks,
# docs/techreport/techreport.tex Sec. "Reproduction", code/scripts/run_e1.sh /
# run_scaling.sh, results/*/run.log), cheapest first. Outputs land in the apps' default
# `results/**/<x>_julia/` directories (the C++ results stay untouched for diffing); one
# log per app under results/logs/; afterwards
#     julia --project=Kirigami/scripts Kirigami/scripts/diff_results.jl --all
# writes results/JULIA_VS_CPP.md.
#
#   nohup bash Kirigami/scripts/rerun_all.sh > results/logs/rerun_all.log 2>&1 &
#   JULIA=/path/to/arm64/julia NSHARDS=12 bash Kirigami/scripts/rerun_all.sh [--only k5,k6] [--dry-run]
#
# Requires an arm64 Julia (bit-exact trig against the C++ build; an x86_64 Julia under
# Rosetta runs everything but tilings/relaxation may round 1 ulp differently, see
# Kirigami/src/core/generators.jl). Estimated wall times are the C++ per-row `secs` sums
# of the archived CSVs (results/*), divided by the shard count where the run is sharded;
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
  echo "WARNING: $JULIA is $ARCH, not arm64 -- the C++ was arm64; tiling/relaxation trig may differ by 1 ulp." >&2
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
run_app() {  # run_app NAME EST "command ..." : one nohup'd process, its own log, wait
  local name=$1 est=$2; shift 2
  selected "$name" || return 0
  if [ ! -f "Kirigami/apps/kill_${name%%[_/]*}.jl" ] && [ ! -f "Kirigami/apps/kill_${name}.jl" ]; then
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
  if [ ! -f "Kirigami/apps/kill_${name%%[_/]*}.jl" ]; then
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
# Cheapest first. Dependencies: k1a -> {k1b, k2c, k8a} (shape cache results/kill/cache,
# shared with the C++ layout); k5 -> k6 (sigma_def, default = corpus) -> {t1, k9, k9b,
# k9c, b4, k8a_recheck} (K6 cache); k9 + k9b -> k9c (their CSVs). Every app defaults to the
# frozen corpus populations, so the order below only matters for the caches.
# ======================================================================================

# --- seconds-scale ---------------------------------------------------------------------
run_app f23      "~1 s"      "$APP/kill_f23.jl"
run_app k2a      "~7 s"      "$APP/kill_k2a.jl --n 292"
run_sharded jitter "~2 s"    --nshards "$APP/kill_jitter.jl"
if selected jitter && [ "$DRY" = 0 ]; then
  merge_csv results/kill/jitter_julia/jitter.csv results/kill/jitter_julia/jitter_[0-9]*.csv
  run_app jitter_analyze "~1 s" "$APP/kill_jitter.jl --analyze results/kill/jitter_julia/jitter.csv"
fi
run_app k3a_recheck "~8 s (C++ 8 s)"   "$APP/kill_k3a_recheck.jl"
run_app k1b      "~14 s (C++ 14 s)"    "$APP/kill_k1b.jl --n 200 --shapes 20"
run_app k1a      "~36 s (C++ 36 s)"    "$APP/kill_k1a.jl --n 200 --samples 10000"
run_app k1c      "~31 s, native CLI"   "$APP/kill_k1c.jl"
run_app k2c      "~47 s (C++ 47 s)"    "$APP/kill_k2c.jl --n 500"
run_app t1       "~7 s off the K6 cache" "$APP/kill_t1.jl"
run_app b3       "~1 min"              "$APP/kill_b3.jl --stage run"

# --- minutes-scale ----------------------------------------------------------------------
run_app k6       "~2 min (C++ 96 s + native)"  "$APP/kill_k6.jl --n 200 --maxf 800"
run_app k5       "~5 min (C++ 307 s)"          "$APP/kill_k5.jl --n 200 --maxf 800"
run_sharded e1   "~6 min / $NSHARDS (C++ 352 s total)" --nshards "$APP/kill_e1.jl --n-random 900"
if selected e1 && [ "$DRY" = 0 ]; then
  merge_csv results/final/e1_julia/e1.csv results/final/e1_julia/e1_[0-9]*.csv
fi
run_sharded k8a  "~11 min / $NSHARDS (C++ 658 s total)" --nshard "$APP/kill_k8a.jl --n 100 --x0 --dual-iters 100000"
if selected k8a && [ "$DRY" = 0 ]; then
  merge_csv results/kill/k8a_julia/k8a.csv results/kill/k8a_julia/k8a_[0-9]*.csv
  run_app k8a_summary "~1 s" "$SCR/summarise_k8a.jl results/kill/k8a_julia"
fi
run_app k8a_recheck "~5 min"           "$APP/kill_k8a_recheck.jl"
run_app k3a      "~24 min (C++ 1410 s)" "$APP/kill_k3a.jl --n 500"
# K7: main + c3 + perdbg sharded (wall set by voronoi_torus_11_n200, ~30 min alone), then
# the three single-process stages
run_sharded k7_main   "~30 min (C++ wall, 12 shards)" --nshard "$APP/kill_k7.jl --stage main"
run_sharded k7_c3     "~30 min"                      --nshard "$APP/kill_k7.jl --stage c3"
run_sharded k7_perdbg "~5 min"                       --nshard "$APP/kill_k7.jl --stage perdbg"
if selected k7_main && [ "$DRY" = 0 ]; then
  merge_csv results/kill/k7_julia/k7_main.csv results/kill/k7_julia/k7_main_[0-9]*.csv
  merge_csv results/kill/k7_julia/k7_c3_all.csv results/kill/k7_julia/k7_c3_[0-9]*.csv
  merge_csv results/kill/k7_julia/k7_periodicity.csv results/kill/k7_julia/k7_periodicity_[0-9]*.csv
fi
run_app k7_bounded "~1 min" "$APP/kill_k7.jl --stage bounded"
run_app k7_nu      "~1 min" "$APP/kill_k7.jl --stage nu"
run_app k7_c4sweep "~1 min" "$APP/kill_k7.jl --stage c4sweep"

# --- hours-scale ------------------------------------------------------------------------
run_sharded b4   "~72 min / $NSHARDS (C++ 4313 s total)" --nshards "$APP/kill_b4.jl --n 200"
if selected b4 && [ "$DRY" = 0 ]; then
  merge_csv results/kill/b4_julia/b4.csv results/kill/b4_julia/shard_[0-9]*.csv
fi
run_sharded k9   "~96 min / $NSHARDS (C++ 5774 s total)" --nshards "$APP/kill_k9.jl --n 200"
run_app k9_aggregate "~1 s" "$APP/kill_k9.jl --aggregate --nshards $NSHARDS"
run_sharded k9b  "~7.7 h / $NSHARDS (C++ 27685 s total)" --nshards "$APP/kill_k9b.jl --n 200"
run_app k9b_aggregate "~1 s" "$APP/kill_k9b.jl --aggregate --nshards $NSHARDS"
run_sharded k9c  "~5.4 h / $NSHARDS (C++ 19318 s total)" --nshards "$APP/kill_k9c.jl --n 200"
run_app k9c_aggregate "~1 s" "$APP/kill_k9c.jl --aggregate --nshards $NSHARDS"
if grep -q "TODO(julia-port)" Kirigami/scripts/run_scaling.sh; then
  echo "[$(stamp)] SKIP scaling: Kirigami/scripts/run_scaling.sh still drives the C++ binary (TODO(julia-port))"
else
  run_app scaling "hours; sequential cells with wall caps (see run_scaling.sh)" \
    "bash Kirigami/scripts/run_scaling.sh results/scaling_julia"
fi
# regime: the C++ ran 4 shards for the constructive arms and 4 for the native cells
REGIME_SHARDS=4
if selected regime; then
  save=$NSHARDS; NSHARDS=$REGIME_SHARDS
  run_sharded regime_pop    "~1 h / 4 (constructive arms)" --nshards "$APP/kill_regime.jl --mode pop"
  run_app     regime_ref    "~10 min"                       "$APP/kill_regime.jl --mode ref"
  NSHARDS=$save
fi
run_app yield    "hours (not yet ported)"  "$APP/kill_yield.jl --mode fresh"
run_app basin    "hours (not yet ported)"  "$APP/kill_basin.jl"

# --- native CLI last (hours; baseline/native/build/tuttekiri_cli must exist) ------------
run_app k2b      "~26 min (C++ 1517 s), native CLI" "$APP/kill_k2b.jl --maxf 160"
# regime native cells (4 shards, 600 s timeout each), then the merge
if selected regime; then
  save=$NSHARDS; NSHARDS=$REGIME_SHARDS
  run_sharded regime_native "hours, native CLI (600 s timeout per cell)" --nshards "$APP/kill_regime.jl --mode native --timeout 600"
  NSHARDS=$save
  run_app     regime_aggregate "~1 s" "$APP/kill_regime.jl --mode aggregate --merge-shards $REGIME_SHARDS"
fi
run_sharded native200 "~65 h / $NSHARDS CPU (C++ 234293 s total, 600 s timeout per cell)" --nshards \
  "$APP/kill_native200.jl --n 200 --maxf 800 --timeout 600 --out results/kill/native200_julia/shards --sigma results/kill/k5_julia/sigma"

# ---- the comparison ------------------------------------------------------------------
if [ "$DRY" = 0 ]; then
  echo "[$(stamp)] diff_results --all"
  $SCR/diff_results.jl --all > "$LOGS/diff_results.log" 2>&1
  echo "[$(stamp)] done -> results/JULIA_VS_CPP.md"
fi
