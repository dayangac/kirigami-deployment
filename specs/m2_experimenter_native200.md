# Mission 2 / WP5 — Experimenter-N: complete the Native200 baseline

Paste of specs/common_preamble.md applies (read it first).

## Why this exists
results/kill/native200/native200.csv holds 573 of 600 (graph, variant) cells of the authors' full native pipeline run on the K9 population: completed 172, timed_out 362 (600 s alarm), crashed 39; 27 cells were never dispatched (the user stopped the run). REPORT.md and the technical report can therefore only say "0 of what finished". This work package finishes the run so the baseline sentence can be printed without that qualifier, or with a quantified one.

## Read first
- results/kill/KILL_REPORT.md §"Native200" (the last section) — how the run was done, the F24 crash mechanism, the false-negative note.
- code/apps/kill_native200.cpp in full (driver; options --n --maxf --out --sigma --cli --work --timeout --shard --nshards --logdir).
- code/apps/kill_common.hpp (population is a deterministic function of graph id; kind = voronoi/delaunay/quad_random by id % 3).
- results/kill/native200/{summary.txt, run_all.log}; the shards/ and logs/ layout.
- STATE.md F24 (upstream bugs: get_holes() segfault from an unchecked prev()->twin(); merge_close_verts fusing hinge duplicates) and F40.
- baseline/native/ (the vendored authors' code and CMakeLists; the CLI is baseline/native/build/tuttekiri_cli).

## Task
1. **Make the driver resumable and cell-selectable.** Add to kill_native200.cpp: `--cells <file>` (a list of `id,variant` lines to run, ignoring the shard split), `--resume <csv>` (skip any cell whose status in that CSV is `completed` or `crashed`), and `--tag <string>` written into a new `run` column so the 600 s run, the 3600 s rerun and the crash-fixed rerun are distinguishable in one merged CSV. Keep the existing per-cell scoring untouched (exact T4.2″ scan, bisection referee, ε = 0.3 certificate, the authors' own native_theta_collisions as a diagnostic column). Rebuild only that target: `cmake --build code/build --target kill_native200`.
2. **Cell lists.** From native200.csv produce results/kill/native200/cells_timed_out.txt (362) and cells_crashed.txt (39); compute cells_missing.txt (the 27 never dispatched) from the full 200 × {native, sigma_mc, sigma_def} grid. Sort the timed-out and missing lists by **|F| ascending** so the cheap cells report first.
3. **Launch the 3600 s rerun DETACHED** — it may take up to 60 h and must survive this agent ending. Write results/kill/native200/run_rerun.sh that runs the timed-out + missing lists at `--timeout 3600 --tag rerun3600`, **6-way parallel** (not 12), writing shard CSVs to results/kill/native200/rerun3600/ and logs to rerun3600/logs/, and a merged progress file rerun3600/progress.txt updated after every cell (count completed / timed_out / crashed so far, and the last id done). Launch with `nohup bash run_rerun.sh > rerun3600/run_all.log 2>&1 &` and record the PID in rerun3600/PID. Verify within a few minutes that cells are actually being produced (watch the smallest-|F| cells complete), then move on; do NOT wait for the run.
4. **Crash-fixed rerun.** In the vendored copy baseline/native, find the unchecked `prev()->twin()` in UnitPattern::get_holes() (F24) and add the minimal null check that makes it return gracefully instead of segfaulting; record the diff in results/kill/native200/crashfix.patch (git diff of baseline/native). Rebuild the CLI into a **separate** binary, e.g. baseline/native/build_fixed/tuttekiri_cli, so the original binary is untouched. Confirm the fix on 3 of the 39 crashed cells manually. Then queue the 39 crashed cells at `--timeout 3600 --tag crashfix3600 --cli baseline/native/build_fixed/tuttekiri_cli` into the SAME detached pipeline (append to run_rerun.sh before launching, run after the timed-out list, or as a second nohup job at 2-way parallelism — your call, but total concurrency must stay ≤ 8 processes; the machine is shared with other agents).
5. **Merge + report tooling.** Write code/apps/native200_merge.cpp (or extend kill_native200 with `--merge`) that combines native200.csv and every shard CSV under rerun3600/ and crashfix3600/ into results/kill/native200/native200_final.csv, one row per cell with the LAST status (a rerun row supersedes the 600 s timed_out row; keep the original row in a `first_status` column). Write results/kill/native200/NATIVE200_FINAL.md from that CSV with:
   - table: status × variant, status × family, for the original run, the 3600 s rerun and the crash-fixed rerun separately and merged;
   - the count of completed cells with exact Θ_max > 0, with bisection referee, and with the ε = 0.3 certificate; the count where the authors' own test says Θ > 0 (false negatives, per the K1c mechanism), listed individually as the existing section does;
   - wall-time distribution of completed cells (median, q90, max) and, for cells that still time out at 3600 s, the |F| above which their pipeline does not finish on this machine;
   - **the single sentence the paper may print about the baseline**, plus the qualifier if any cells remain unfinished.
   The report must be regenerable: NATIVE200_FINAL.md is produced by a script (C++ writer or a ≤ 100-line Python that only reads the CSV) so the orchestrator can rerun it as the detached job progresses. Generate a first version now, labelled "PARTIAL — rerun in progress, N of 389 cells done".
6. **Do not** rerun the 172 completed cells, and do not touch results/kill/k9c or any other results directory.

## Rules
- Concurrency ≤ 8 native processes in total. Check `uptime` before launching.
- Every number in NATIVE200_FINAL.md comes from native200_final.csv.
- Do not commit. Do not edit REPORT.md/IDEA.md/STATE.md.
- If the CLI's behaviour on the crash-fixed build differs from the original on a *completed* cell (compare 3 completed cells' Θ output between the two binaries), report that; it decides whether the crash-fixed rows can be pooled with the rest.

## Reply
≤15 lines: what was changed in the driver; the PID and the command to check progress; how many cells had completed when you left; the crash-fix diff in one line and whether it changes completed-cell output; the path of the report generator and how to rerun it.
