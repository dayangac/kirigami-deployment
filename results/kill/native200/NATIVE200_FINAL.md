# Native200 -- final status of the authors' full native pipeline on the K9 population

**PARTIAL -- rerun in progress, 4 of 389 cells done** (3600 s rerun of the timed-out and never-dispatched cells), plus 26 of 39 crash-fixed cells. Regenerate with `code/build/native200_merge`.

Every number below is computed by `Kirigami/apps/native200_merge.jl` from `results/kill/native200/native200_final.csv`, which is itself merged from `native200.csv` (the 600 s run), `rerun3600/shard_*.csv` (the 3600 s rerun of the timed-out and never-dispatched cells) and `crashfix3600/shard_*.csv` (the 39 crashed cells against the F24-patched CLI). One row per (graph id, variant); the 600 s status of a superseded cell is kept in `first_status`.

Grid: 600 cells (200 graphs x {native, sigma_mc, sigma_def}). Rows by run: crashfix3600 26; never dispatched 26; orig (600 s) 544; rerun3600 4; 

## 1. Status by variant

**ALL**

| | native | sigma_mc | sigma_def | total |
|---|--:|--:|--:|--:|
| completed | 58 | 64 | 73 | 195 |
| crashed | 1 | 4 | 10 | 15 |
| missing | 11 | 7 | 8 | 26 |
| timed_out | 130 | 125 | 109 | 364 |
| **total** | **200** | **200** | **200** | **600** |

**crashfix3600**

| | native | sigma_mc | sigma_def | total |
|---|--:|--:|--:|--:|
| completed | 1 | 0 | 21 | 22 |
| crashed | 0 | 0 | 2 | 2 |
| timed_out | 0 | 2 | 0 | 2 |
| **total** | **1** | **2** | **23** | **26** |

**never dispatched**

| | native | sigma_mc | sigma_def | total |
|---|--:|--:|--:|--:|
| missing | 11 | 7 | 8 | 26 |
| **total** | **11** | **7** | **8** | **26** |

**orig (600 s)**

| | native | sigma_mc | sigma_def | total |
|---|--:|--:|--:|--:|
| completed | 57 | 64 | 51 | 172 |
| crashed | 1 | 4 | 8 | 13 |
| timed_out | 129 | 121 | 109 | 359 |
| **total** | **187** | **189** | **168** | **544** |

**rerun3600**

| | native | sigma_mc | sigma_def | total |
|---|--:|--:|--:|--:|
| completed | 0 | 0 | 1 | 1 |
| timed_out | 1 | 2 | 0 | 3 |
| **total** | **1** | **2** | **1** | **4** |

## 2. Status by family

**ALL**

| | delaunay | voronoi | quad_random | total |
|---|--:|--:|--:|--:|
| completed | 120 | 46 | 29 | 195 |
| crashed | 5 | 1 | 9 | 15 |
| missing | 9 | 7 | 10 | 26 |
| timed_out | 67 | 147 | 150 | 364 |
| **total** | **201** | **201** | **198** | **600** |

**crashfix3600**

| | delaunay | voronoi | quad_random | total |
|---|--:|--:|--:|--:|
| completed | 2 | 4 | 16 | 22 |
| crashed | 0 | 0 | 2 | 2 |
| timed_out | 1 | 0 | 1 | 2 |
| **total** | **3** | **4** | **19** | **26** |

**never dispatched**

| | delaunay | voronoi | quad_random | total |
|---|--:|--:|--:|--:|
| missing | 9 | 7 | 10 | 26 |
| **total** | **9** | **7** | **10** | **26** |

**orig (600 s)**

| | delaunay | voronoi | quad_random | total |
|---|--:|--:|--:|--:|
| completed | 118 | 41 | 13 | 172 |
| crashed | 5 | 1 | 7 | 13 |
| timed_out | 64 | 146 | 149 | 359 |
| **total** | **187** | **188** | **169** | **544** |

**rerun3600**

| | delaunay | voronoi | quad_random | total |
|---|--:|--:|--:|--:|
| completed | 0 | 1 | 0 | 1 |
| timed_out | 2 | 1 | 0 | 3 |
| **total** | **2** | **2** | **0** | **4** |

## 3. What the completed cells say

| quantity | count | of completed |
|---|--:|--:|
| completed cells (any run) | 195 | 195 |
| exact `Theta_max > 0` (T4.2" scan on their dumped embedding) | 0 | 195 |
| bisection referee `Theta > 0` | 0 | 195 |
| validity certificate at eps = 0.3 (POS & NOOVERLAP & NOROOT) | 0 | 195 |
| the authors' own FK collision test says `Theta > 0` | 27 | 195 |

The cells where their own test reports `Theta > 0` while our exact scan on the SAME dumped embedding gives 0 (the F24/K1c false-negative mechanism):

| id | kind | variant | F | dim_null | native_theta | our_exact | our_bisect | cert_pos | run |
|---|---|---|--:|--:|--:|--:|--:|--:|---|
| 19 | delaunay | sigma_def | 102 | 51 | 2.462 | 0.000 | 0.000 | 0 | crashfix3600 |
| 26 | quad_random | sigma_def | 255 | 206 | 0.777 | 0.000 | 0.000 | 0 | crashfix3600 |
| 36 | voronoi | sigma_def | 131 | 235 | 0.673 | 0.000 | 0.000 | 0 | orig (600 s) |
| 38 | quad_random | sigma_def | 187 | 149 | 1.096 | 0.000 | 0.000 | 0 | crashfix3600 |
| 42 | voronoi | sigma_def | 125 | 230 | 1.091 | 0.000 | 0.000 | 0 | crashfix3600 |
| 44 | quad_random | sigma_def | 319 | 263 | 0.712 | 0.000 | 0.000 | 0 | crashfix3600 |
| 48 | voronoi | sigma_def | 230 | 418 | 0.732 | 0.000 | 0.000 | 0 | crashfix3600 |
| 54 | voronoi | native | 139 | 158 | 0.888 | 0.000 | 0.000 | 0 | crashfix3600 |
| 59 | quad_random | sigma_def | 329 | 268 | 0.272 | 0.000 | 0.000 | 0 | crashfix3600 |
| 65 | quad_random | sigma_def | 230 | 186 | 0.872 | 0.000 | 0.000 | 0 | crashfix3600 |
| 68 | quad_random | sigma_def | 269 | 222 | 0.828 | 0.000 | 0.000 | 0 | crashfix3600 |
| 71 | quad_random | sigma_def | 257 | 208 | 0.983 | 0.000 | 0.000 | 0 | crashfix3600 |
| 72 | voronoi | sigma_def | 200 | 362 | 0.228 | 0.000 | 0.000 | 0 | crashfix3600 |
| 73 | delaunay | sigma_def | 146 | 74 | 1.555 | 0.000 | 0.000 | 0 | crashfix3600 |
| 102 | voronoi | sigma_def | 336 | 600 | 0.351 | 0.000 | 0.000 | 0 | orig (600 s) |
| 128 | quad_random | sigma_def | 252 | 203 | 0.884 | 0.000 | 0.000 | 0 | crashfix3600 |
| 146 | quad_random | sigma_def | 187 | 151 | 0.799 | 0.000 | 0.000 | 0 | crashfix3600 |
| 152 | quad_random | sigma_def | 290 | 243 | 0.915 | 0.000 | 0.000 | 0 | crashfix3600 |
| 155 | quad_random | sigma_def | 298 | 234 | 0.975 | 0.000 | 0.000 | 0 | crashfix3600 |
| 156 | voronoi | sigma_def | 308 | 566 | 0.570 | 0.000 | 0.000 | 0 | orig (600 s) |
| 160 | delaunay | native | 292 | 59 | 0.088 | 0.000 | 0.000 | 0 | orig (600 s) |
| 168 | voronoi | sigma_def | 138 | 238 | 1.025 | 0.000 | 0.000 | 0 | orig (600 s) |
| 173 | quad_random | sigma_def | 213 | 177 | 0.573 | 0.000 | 0.000 | 0 | crashfix3600 |
| 176 | quad_random | sigma_def | 188 | 157 | 1.016 | 0.000 | 0.000 | 0 | crashfix3600 |
| 179 | quad_random | sigma_def | 188 | 152 | 0.345 | 0.000 | 0.000 | 0 | crashfix3600 |
| 185 | quad_random | sigma_def | 274 | 218 | 0.719 | 0.000 | 0.000 | 0 | crashfix3600 |
| 197 | quad_random | sigma_def | 194 | 157 | 0.970 | 0.000 | 0.000 | 0 | crashfix3600 |

## 4. Cost

| quantity | value |
|---|--:|
| completed wall time, median (s) | 13.480 |
| completed wall time, q90 (s) | 186.612 |
| completed wall time, max (s) | 570.706 |
| completed \|F\|, median | 202.000 |
| completed \|F\|, max | 760.000 |
| cells still timing out at 3600 s | 5 |
| smallest \|F\| that still times out at 3600 s | 128.000 |
| median \|F\| of the cells that still time out | 131.000 |

On this machine their pipeline does not finish within 3600 s above roughly |F| = 128.000 faces (the smallest patch that still times out); the largest patch that does finish has |F| = 760.000.

## 5. The sentence the paper may print

> On the 200-graph K9 population, the authors' full native pipeline (their colouring, their Eq. (6), their Eq. (9) at the published web-UI defaults, their forward-kinematics collision test) reaches an exact `Theta_max > 0` on **0 of 195** designs it completes, and passes the eps = 0.3 validity certificate on **0 of 195**; their own collision test reports `Theta > 0` on 27 of them, every one a false negative of the F24/K1c merge_close_verts mechanism.

Qualifier that must accompany it: **405 of 600 cells (67.5 %)** never produce a design at all -- they time out or crash -- so the denominator is the completed subset, not the population.
