# F32 — the NOROOT false negatives: diagnosis, fix, and re-measurement

Author: subagent (cert-diagnosis). Everything below is measured on this machine with the
Julia code in `Kirigami/`; no number is quoted from memory.

## 1. Symptom (as reported in STATE.md F32)

`validity_certificate(..., eps = 0.006)` returned `noroot = false` on 1622/1859
split-bearing jitter rows whose exact `Theta_max` (the T4.2'' interval scan
`exact_theta_max_overlap`) was 1.0–2.4 rad. Re-measured on the pre-fix CSVs kept in
`results/kill/jitter/pre_f32fix/`, the same effect with the metric restated as
"split-bearing rows with `Theta_max >= 1` that the certificate does not certify":

| population (6416 solved rows, pre-fix) | count |
|---|---|
| certified at `eps = 0.006` | 3445 (53.7%) |
| split-bearing rows with `Theta_max >= 1` | 1161 |
| …of those, NOT certified | 1010 (87.0%) |
| rows with exact `Theta_max >= eps`, not certified | 1619 / 5064 (32.0%) |
| rows certified but with `Theta_max < eps` (unsound) | 0 |

So the certificate was sound and badly incomplete.

## 2. Root cause: a missing predicate, not a tolerance

**The certificate's NOROOT tested a strictly larger root set than the exact scan's
contact set `C(X)`. It omitted the two interval (projection) inequalities of T4.1b.**

A root of the orientation harmonic `h_o,pi(theta)` only says the vertex `w` is collinear
with the **infinite line** through the edge `(a,b)`. A contact — and hence a transition of
the overlap status, which is all Proposition T5.2b' needs — additionally requires `w` to
lie on the **segment**:

```
        0  <=  <w - a, b - a>  <=  |b - a|^2 .
```

`contact_angles` in `Kirigami/src/method/contact.jl` (the exact scan) applies exactly that test, with
`s = D.eval(th)`, `l2 = L2.eval(th)`, `tol = 1e-12 * lscale`.
`validity_certificate`'s NOROOT scan (same file) did not: it called
`harmonic_roots_deflated(det, 0, eps)` and set `noroot = false` on any root.

This was not an implementation slip. `derivations/core.md` T5.2b.0 states it deliberately:

> `(iii)` is stricter than the contact condition — it drops the two interval tests `E2`
> and the overlap selection — so `R(eps)` is smaller than the exact region.

The 87% rejection rate is the measured price of that drop.

### Evidence: five rows, every root listed

`Kirigami/apps/dbg_cert.jl` (new) re-runs one jitter row through the identical `measure()`
pipeline and prints, for every root the certificate finds in `(0, eps)`, the ordered pair,
the root, the deflation class, `h(0)`, the projection parameter `s/|e|^2`, whether the
vertex is inside the segment, and whether the root is a graze (`|h'(theta)|` below
`1e-9 * scale`).

```
julia --project=Kirigami Kirigami/apps/dbg_cert.jl <tiling> <amp_idx> <seed> [eps]   # diagnostic of the pre-F32 pass; the app is not kept
```

| row (tiling, amp_idx, seed, sigma_mc) | exact `Theta_max` | roots in `(0,eps)` | inside segment | outside | graze |
|---|---|---|---|---|---|
| `hexagons_auto` 0 0 (a = 0.005) | 1.027137 | 32 | **0** | 32 | 0 |
| `hexagons_auto` 0 1 (a = 0.005) | 1.023071 | 33 | **0** | 33 | 0 |
| `snub_square_33434` -1 0 (a = 0, the authored tiling) | 1.646136 | 2 | **0** | 2 | 0 |
| `tiling_3_4_3_12` -1 0 (a = 0, the authored tiling) | 2.380636 | 4 | **0** | 4 | 0 |
| `truncated_square_488` 0 0 (a = 0.005) | 2.332072 | 97 | **0** | 97 | 0 |

All 168 reported roots have the vertex off the segment. Representative lines:

```
truncated_square_488 ai=0 seed=0  Theta_max = 2.332072
  root th=3.551501e-03 faces(5->0) v=21 edge(3,0) class=c1 h(0)=8.537e-03  s/|e|^2= 3.4234 in_seg=0 graze=0
  root th=3.973755e-03 faces(3->0) v=14 edge(0,1) class=c2 h(0)=0.000e+00  s/|e|^2=-0.0039 in_seg=0 graze=0
tiling_3_4_3_12 ai=-1 seed=0     Theta_max = 2.380636
  root th=3.388164e-03 faces(3->14) v=14 edge(49,50) class=c1 h(0)=2.815e-03 s/|e|^2= 6.7767 in_seg=0 graze=0
```

The projection parameter runs to `6.78` (nearly seven edge lengths past the far endpoint)
on 3.4.3.12 and to `-3.0`/`+4.0` on hexagons. The class-2 roots sit at
`s/|e|^2 ~ -0.004`, just past the near endpoint.

### The four candidate causes, each tested

| candidate | verdict | evidence |
|---|---|---|
| (a) near-zero roots from permanent incidences, `h(0) ~ 1e-12 != 0` not struck | **NO** | `n_identically_zero = 0` on all five rows; the class-2 roots have `h(0)` exactly `0.000e+00` and are correctly deflated (their reported root is the *other* root `2 atan(-B/A)`, not the flat state). The class-1 roots have `h(0)` of order `1e-3`, i.e. `~1e-3 * scale` — nowhere near the identity tolerance. |
| (b) graze (double) roots admitted as crossings | **NO** | 0 grazes among the 168 roots; every root has `\|h'\|` far above `1e-9 * scale`. |
| (c) certificate and exact scan use different predicate sets | **YES — this is the cause** | The exact scan filters by the interval inequalities, the certificate did not. 168/168 of the reported roots fail that filter; `contact_angles(c, B, pairs, 0, eps)` is empty on all five rows while the certificate reported 2–97 roots. |
| (d) `eps` handling / `tan(theta/2)` chart | **NO** | `harmonic_roots` uses the amplitude/phase form (`atan2`+`acos`), not the `tau` chart (check.md D10); the roots reported are genuine roots of the harmonic and are reproduced independently by `dbg_cert`. |

## 3. The fix

`Kirigami/src/method/contact.jl`, in `validity_certificate`'s NOROOT scan: build the edge's
squared-length harmonic `L2` and the projection harmonic `D` exactly as `contact_angles`
does, and keep only the roots that satisfy `-tol <= D.eval(th) <= L2.eval(th) + tol` with
`tol = 1e-12 * (|L2.p| + L2.amp())`. The test is **inclusive at the tolerance**, so a
borderline root stays admissible and the certificate stays on the conservative side.
Roots removed by the filter are counted in a new field
`ValidityCertificate.n_roots_inadmissible` (reporting only, never a violation).

Nine lines of substance; no other file's logic changes.

**Why this is sound.** After the change, NOROOT(eps) is exactly the statement
`C(X) ∩ (0, eps) = ∅` for `C(X)` the contact-angle set of T4.1b / Corollary T4.2'. The
connectedness argument of Proposition T5.2b' needs precisely that: the overlap status can
change only at an angle of `C(X)` (T4.2'), so if `C(X)` misses `(0, eps)` the status is
constant there, and `NOOVERLAP(eps/2)` fixes the constant. Dropping the interval tests was
never needed by the proof; it only shrank `R(eps)`.

## 4. Proposed amendment to `derivations/core.md` T5 (NOT applied — core.md is untouched)

In **T5.2b.0(iii)**, replace

> `(iii) NOROOT(t) : for every candidate pair pi in the C-list, the orientation harmonic
> h_o,pi has no root in (0, eps)`

by

> `(iii) NOROOT(t) : for every candidate pair pi = (w,(a,b)) in the C-list, the
> orientation harmonic h_o,pi has no ADMISSIBLE root in (0, eps), where a root theta is
> admissible iff it also satisfies the two interval inequalities of T4.1b,`
> ```
>       0  <=  <w(theta) - a(theta), b(theta) - a(theta)>  <=  |b(theta) - a(theta)|^2 ,
> ```
> `both sides being harmonics in theta by T3.4. Equivalently: C(X) ∩ (0, eps) = ∅.`

and delete the sentence

> "`(iii)` is stricter than the contact condition — it drops the two interval tests `E2`
> and the overlap selection — so `R(eps)` is smaller than the exact region."

replacing it with: "`(iii)` drops only the overlap *selection* of T4.2''; it keeps the
interval tests, so `R(eps)` is smaller than the exact region only through
`NOOVERLAP(theta_1)` being tested at one angle."

**One caveat a Checker must adjudicate.** Sub-lemma T5.2b'' Case A (round 5) handles the
case where the Lemma-T4.2 witness at the limit point is a permanent incidence, by
producing a *substitute* candidate pair — "far endpoint of the shorter far-side edge on
the other far-side edge" — whose harmonic is `± L L' sin(beta_e - theta)` and vanishes at
`beta_e`. The printed proof asserts that this pair is in the candidate list and that its
harmonic vanishes; it does **not** assert that `beta_e` is an *admissible* root of that
pair, i.e. that the far endpoint lies on the segment at `theta = beta_e`. Under the
amended (iii) that step now needs the extra sentence. Geometrically it holds: at
`theta = beta_e` the two far-side edges are collinear and meet at the far-side gap
closing, so the endpoint is an endpoint of the overlap of the two collinear segments, and
the interval predicate holds with equality at worst — which the inclusive tolerance
admits. It is stated here as the one gap the amendment opens, not as a proved step.
`derivation_tests.jl`' R5-a…R5-d block (which checks (T5.2b''-1b) on the corpus) passes
unchanged after the fix, but it does not test this admissibility point.

## 5. Verification

**Regression test** (`Kirigami/test/test_method.jl`, "F32: NOROOT counts only roots with the
vertex ON the edge segment"): builds 4.8.8, snub square and 3.4.3.12 at their unjittered
positions, keeps those with exact `Theta_max >= 1`, and requires `POS`, `NOOVERLAP`,
`NOROOT`, `n_roots_deflated == 0`, `contact_angles(...,0,eps)` empty, and
`n_roots_inadmissible > 0` summed over the cases.

* Without the fix (`git stash` of `contact.jl` only, everything else identical):
  `12 assertions, 7 passed, 5 failed` — `cert.noroot` false, `n_roots_deflated != 0`,
  `inadmissible_total == 0`.
* With the fix: passes. Reported off-segment roots removed: 2 (snub square,
  `Theta_max = 1.39352`), 4 (3.4.3.12, `Theta_max = 2.38064`).

Two pre-existing certificate tests had to be updated because their own reference rescans
reproduced the *unfiltered* root set; the same interval filter was added to the rescan, and
their `eps` was widened from 1.2 to 3.0 with two more patterns (snub square, 3.4.3.12)
added, because on the regular patterns almost no root is admissible below 3 rad. These are
changes to the reference computation in the test, not weakenings of the assertions: the
"first_root is the minimum over all candidates" test still finds a discriminating case
(min root 1.0472 vs first-seen 2.0944).

**Suites** (`julia --project=Kirigami -e 'using Pkg; Pkg.test()'`):

| suite | result |
|---|---|
| `Pkg.test()` (whole suite) | 187 test sets / 176,692 assertions (23 marked broken), SUCCESS |
| `derivation_tests.jl` | 38 test sets / 150,191 assertions, SUCCESS |

## 6. Re-measurement of A3 (`Kirigami/apps/kill_jitter.jl`, 12 shards, ~2 s)

Same 6416 solved rows, same geometry, only the certificate changed.

| metric (eps = 0.006) | before | after |
|---|---|---|
| certified | 3445 (53.7%) | **5064 (78.9%)** |
| certified, split-bearing rows only | 237 / 3208 (7.4%) | **1856 / 3208 (57.9%)** |
| `NOROOT` holds | 3580 | **6197** |
| split-bearing rows with `Theta_max >= 1` NOT certified | 1010 / 1161 (**87.0%**) | **0 / 1161 (0.0%)** |
| rows with `Theta_max >= eps` NOT certified | 1619 / 5064 (32.0%) | **0 / 5064 (0.0%)** |
| certified but `Theta_max < eps` (unsound) | 0 | **0** |

**The certificate is now exact on this corpus**: the 5064 certified rows are precisely the
5064 rows with exact `Theta_max >= eps`, in both directions, with no unsound row. The
binding clause is now `NOOVERLAP` (5065 of 6416), as it should be — `NOROOT` holds on 6197.

The A3 verdict changes as a result: `a*(certified)` now coincides with
`a*(Theta_max > 0)` on every tiling and rule (e.g. hexagons mc 0.2722 vs 0.2773; snub
square mc 0.1632 vs 0.1632; 3.4.3.12 mc 0.5217 vs 0.5217; 4.8.8 mc 0.3196 vs 0.3196), so
the transition is a property of the **geometry**, not of the certificate. The sentence
"the certificate, not the geometry, is knife-edge" is withdrawn. Predictor AUCs are
essentially unchanged (`badq_ini` 0.9872 all rows / 0.9760 split-only for
`Theta_max > 0`), so the A3 PASS/FAIL logic on the predictor bar is unaffected.

New `summary.txt` and `ladder.csv` are regenerated in `results/kill/jitter/`; the pre-fix
CSVs and summary are preserved under `results/kill/jitter/pre_f32fix/`.

## 7. Impact on earlier facts

**K2a — the numbers change, the claims survive and get stronger.** Re-ran
`julia --project=Kirigami Kirigami/apps/kill_k2a.jl` (6.5 s, 187 configurations):

| K2a line | before (KILL_REPORT L408–418) | after |
|---|---|---|
| certificate holds, `eps = 0.006` | 65 / 187 | **173 / 187** |
| `NOROOT` alone | (not separately quoted) | **187 / 187** |
| **T5.2b': certificate ⟹ bisection `Theta_max >= eps`** | 65 / 65, 0 violations | **173 / 173, 0 violations** |
| configurations that actually have `Theta_max >= eps` | 173 / 187 | 173 / 187 (unchanged — geometry) |
| spurious roots the `tau = 0` deflation removes | 13 630 | 13 630 (unchanged) |

The soundness claim is unaffected in direction and is now tested on 173 configurations
instead of 65. The sentence "the certificate is a strict **inner** approximation — 65
certified against 173 actually valid" must be **replaced**: on K2a's corpus the
certificate is now *exactly* the valid set, 173 = 173. Any text saying NOROOT's dropped
interval tests are what makes it inner is now wrong.

**K5 and K6 — verdicts unaffected at the time of this fix.** Both FAIL because the exact scan finds a `0+`
split-duplicate collision, i.e. exact `Theta_max = 0`; the binding certificate clause
there is `NOOVERLAP(eps/2)` (and `POS`), not `NOROOT`. A weaker `NOROOT` cannot turn a
`Theta_max = 0` design into a certified one. Re-run confirmation is recorded in §8.

**Anything quoting "certified count" as a design-space size** must be re-derived: the
counts above are all larger. F32's standing instruction ("until fixed, no certified count
may be quoted") is discharged for `NOROOT`; the remaining conservatism of the certificate
is `NOOVERLAP` at the single angle `eps/2`, which on these corpora is exact.

## 8. K5 / K6 re-runs after the fix (measured, not argued)

**K5** (`julia --project=Kirigami Kirigami/apps/kill_k5.jl`, full 200 graphs, both sigma rules, `eps = 0.3`) — identical to the
recorded numbers:

| | `sigma_mc` | `sigma_def` |
|---|---|---|
| valid `X0` by the certificate | **0 / 200** | **0 / 200** |
| POS | 58 | 186 |
| NOOVERLAP(`eps/2`) | 0 | 0 |
| NOROOT(`eps`) | 0 | 0 |
| exact `Theta_max == 0` | 200 / 200 | 200 / 200 |
| certified `Theta_max > 0` | 0 | 0 |

Verdict on the final rule: **FAIL**, unchanged. NOROOT is still 0/200 after the fix,
because the `0+` split-duplicate collision is a *genuine* contact: the vertex really is on
the segment there, so the admissibility filter does not remove it. That is exactly the
discrimination the fix buys — it removes off-segment collinearities and keeps real
contacts.

**K6** (`Kirigami/apps/kill_k6.jl --n 20 --no-native`, a smoke re-run under the same code): every clause
still fails on every graph, `eps_max = 0` and exact `Theta_max = 0` throughout, so the
"certified `Theta_max > 0` on 0/400" verdict is unaffected. Counts from the completed
smoke run are in §8b.

**Summary of impact on earlier facts.**

| earlier fact | status after the fix |
|---|---|
| K2a "certificate implies `Theta_max >= eps` on **65/65**, 0 violations" | **superseded, in the same direction**: now **173/173**, 0 violations |
| K2a "173/187 configurations actually have `Theta_max >= eps`" | **unchanged** (a property of the geometry) |
| K2a "the certificate is a strict inner approximation — 65 certified against 173 valid" | **WRONG after the fix**: 173 certified against 173 valid, i.e. exact on this corpus |
| K5 FAIL (0/200 both sigma) | **unchanged** |
| K6 FAIL (0/400 certified `Theta_max > 0`) | **unchanged by this fix** (the count moved later for other reasons: 2/400 on the primary variant after the referee's collision predicate was corrected, `results/core_validation/referee_fix.md`; 17/400 deployable and 15 certified at `eps_max >= 0.1` rad under the K6 headline rule, best of the four repair variants, `results/kill/KILL_REPORT.md` §K6 — still a FAIL against the 20 % bar) |
| A3 "certificate rejects 1622/1859 rows with `Theta_max >= 1`; the certificate, not the geometry, is knife-edge" | **WITHDRAWN**: 0 rejections after the fix; the transition is geometric |
| F32's ban on quoting certified counts | **discharged for NOROOT**; certified counts must all be re-measured, and every one of them rises |

### 8b. K6 smoke re-run, completed

`julia --project=Kirigami Kirigami/apps/kill_k6.jl --n 20 --no-native` (40 rows = 20 graphs x 2 sigma):

| clause / quantity | count |
|---|---|
| certificate at `eps = 0.3` | **0 / 40** |
| POS | 1 / 40 |
| NOOVERLAP(`eps/2`) | 0 / 40 |
| NOROOT(`eps`) | 0 / 40 |
| `eps_max > 0` | **0 / 40** |
| exact `Theta_max == 0` | **40 / 40** |

Identical to the recorded K6 outcome. The K6 FAIL verdict and its "certified
`Theta_max > 0` on 0/400" are unaffected by this fix, as expected: the binding failure is
a real `0+` contact that the admissibility filter correctly keeps.
