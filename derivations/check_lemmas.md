# check_lemmas.md — adversarial verification of `derivations/lemmas.md` (Checker-L)

**Order of work, as required by the spec, and as actually performed.**

1. I read `specs/m2_checker_lemmas.md` and `specs/common_preamble.md`.
2. **Pass 1 was blind.** The only mathematics I read before writing the Pass-1 verdicts below
   was `derivations/lemmas_statements.md` (statements + standing hypotheses only),
   `derivations/core.md` §0, T1, T3, T5.1 + the deflation block (T5.2b.1/T5.2b.2/Lemma T5.1e/
   the three-class atom list), T5.3, `derivations/check.md` entries D5, R2.5, R3.1a and Round 6,
   and the code (`deploy_basis.hpp`, `contact.hpp`, `zero_plus.hpp`, `periodic_jacobian.hpp`,
   `apps/kill_b3.cpp::achievable/shape_point`, `tests/derivation_tests.cpp`).
   **`derivations/lemmas.md` was NOT opened**, by me or by any tool, until this Pass-1 section
   had been written to disk in full. Every algebraic identity in Pass 1 was re-derived from
   (T1.5)/(T3.2) by hand in the code's sign convention (0.7).
3. Only then Pass 2 (comparison against `lemmas.md`), then the tests.

Verdict scale: PROVED-INDEPENDENTLY / PLAUSIBLE-UNPROVED / COUNTEREXAMPLE / STATEMENT-ILL-POSED.

---

# PASS 1 — blind

Notation as in core.md §0: `c = cos(θ/2)`, `s = sin(θ/2)`, `J` = rotation by `+π/2`,
`C_{(v,f)} = x_v`, `S_{(v,f)} = J(2u_f − σ_f x_v)` (T1.5). Two identities used repeatedly:

```
        det(J a, b) = −⟨a, b⟩ ,        det(a, J b) = ⟨b, a⟩ ,        ⟨a, J b⟩ = −det(a,b).   (K1)
```
(Direct: `Ja = (−a_y, a_x)`.)

## L1.1 — verdict (i) PROVED-INDEPENDENTLY · (ii) PROVED-INDEPENDENTLY for the implication, **GAP in "exactly four types"** · (iii) PROVED-INDEPENDENTLY (and H-L1 is not needed for what (iii) states) · (iv) NOT-CHECKED (measurement)

**(i).** From (T3.2), `p + q = det(C_ab, C_aw)`. By (T1.5) `C_{(v,f)} = x_v` regardless of `f`, so
`C_ab = x_{v_b} − x_{v_a}` and `C_aw = x_{v_w} − x_{v_a}`, giving (L1.6) verbatim. Also
`h_π(0) = p + q·1 + r·0 = p + q`, consistently. `X(t) = X₀ + Φ t` is affine and `det` is bilinear,
so `h_π(0;t)` is a polynomial of degree `≤ 2` in `t`. Dependence only on `(v_a,v_b,v_w)` is
immediate because `C` forgets the face index. **PROVED.**

**(ii) the implication.** `v_w = v_a ⟹ C_aw = 0 ⟹ det = 0`; `v_w = v_b ⟹ C_aw = C_ab ⟹ det = 0`.
Holds for every `X ∈ R^{N×2}`, no hypothesis at all. **PROVED.**

**(ii) the classification "exactly four combinatorial types" — GAP (my id P1-L1-a).**
A coincident-copy pair needs `v := v_w` to be a corner of both `f` and `g`, `f ≠ g`. The listed
types split on *how f and g meet at v*: (a)(b) a shared hinge edge `e ∋ v`, (c) a shared split
edge `e ∋ v`, (d) "share the vertex `v_w` but no edge". The dichotomy is exhaustive only if the
alternative to (a)–(c) is "`f`,`g` share **no edge containing `v_w`**". As printed, (d) says they
share *no edge at all*, which leaves uncovered the configuration
> `f` and `g` are adjacent across an edge `e`, and additionally meet at a vertex `v_w ∉ e`.
That configuration is not excluded by H1–H6 (it needs a non-convex face wrapping onto a pinch
vertex; it is rare but not forbidden by any stated hypothesis), and in it the pair is a genuine
coincident-copy pair belonging to none of (a)–(d) as worded. The repair is one word: (d) should
read "share the vertex `v_w` but no **edge containing `v_w`**", after which the four types are
exhaustive by inspection (the shared edge at `v` is hinge or split; if hinge, `v` is `src(e)` or
`dst(e)`; if split, both endpoints are duplicated). I could not construct a corpus witness, so I
grade this **gap/wording, not wrong**. My test T-L1-c below counts the configuration on the corpus.

Two further remarks I checked and found sound: type (a) is not merely a "permanent incidence" —
under the cut of §0.2 the hinge is *kept* at `src(e)`, so `w` and `a` are literally the **same**
`M′`-vertex, and `h_π ≡ 0` follows from `Y_w ≡ Y_a`, not merely from `h(0) = 0`. And a pair may be
coincident-copy on **both** endpoints simultaneously (type (c) with `v_w = a_e` and the tested edge
being `e` itself); the four types are not claimed to be disjoint, and they are not.

**(iii).** `h_π(0;t)` is a polynomial of degree `≤ 2` on `R^{2k}` (by (i)). Either it is the zero
polynomial, or its zero set is a proper algebraic subset, which for a non-zero real polynomial has
Lebesgue measure zero (standard: the zero set of a non-zero polynomial is a null set). The stated
dichotomy is therefore **unconditionally true** and needs no H-L1; H-L1 is precisely the extra
assertion that the first branch never occurs. The label "(converse, under H-L1)" is therefore a
mild mis-labelling — the *statement* holds without H-L1, only the *usefulness* (that the accidental
set is the whole story) needs it. **PROVED** as stated. **(wording, my id P1-L1-b.)**

**(iv).** A measurement. I did not re-run it and do not certify it. Internal arithmetic is
consistent: `135 410 + 38 944 = 174 354`; `174 354 + 22 251 + 1 116 = 197 721`;
`197 721 / 30 119 368 = 0.66 %`; `1 116 / 30 119 368 = 3.7e−5`. Note that L1.4(iii) quotes
**0.58 %** for the H-L1 failure rate while (iv) headlines **0.66 %**; these are different
quantities (`174 354/30 119 368 = 0.579 %` are the *structural* vanishings, which is the right
number for "H-L1 is false", while `0.66 %` includes the accidental and the `dim_null = 0` rows).
Consistent, but the two numbers should be labelled where they appear. **(wording, P1-L1-c.)**

## L1.2 — verdict PROVED-INDEPENDENTLY (all of (L1.8), the trichotomy, L1.2a, L1.2b); L1.2c PLAUSIBLE-UNPROVED

**The separation identity.** `C_w = x_v = C_{a*}` since both are copies of the same `v`, so
`Y_w − Y_{a*} = c·0 + s·(S_w − S_{a*}) = sin(θ/2)·dS`, exactly, for all `θ`. This is the same
identity `zero_plus.hpp` states as `dC_e = 0` and it is here proved for hinge- and vertex-only
incidences too, not only for split edges. **PROVED.**

**(L1.8), case `a* = a` (`v_w = v_a`).** `C_aw = 0`, `S_aw = dS`. (T3.2) gives
`p = ½ det(S_ab, dS)`, `q = −p`, `r = ½ det(C_ab, dS) = ½ det(d, dS)`.
Within one face (T3.3) `S_ab = −σ_f J C_ab = −σ_f J d`, so by (K1)
`p = −½σ_f det(Jd, dS) = +½ σ_f ⟨d, dS⟩`. **Exactly (L1.8).**

**(L1.8), case `a* = b` (`v_w = v_b`).** Now `C_aw = C_ab = d` and `S_aw = dS + S_ab`. Then
`p = ½[det(d,d) + det(S_ab, dS + S_ab)] = ½ det(S_ab, dS)` — the same expression — `q = −p`, and
`r = ½[det(d, dS + S_ab) + det(S_ab, d)] = ½ det(d, dS)`, the cross terms cancelling. **Same
formulas, with `d` still `x_{v_b} − x_{v_a}` and not re-oriented towards `v`.** This is the one
place the statement could have hidden a sign error and it does not. **PROVED.**

Hence `C = p + q = 0`, `A = p − q = 2p = σ_f ⟨d, dS⟩`, `B = 2r = det(d, dS)`. **PROVED.**

**The trichotomy.** `dS = 0 ⟹ p = q = r = 0 ⟹ h ≡ 0`. Given `dS ≠ 0` and `d ≠ 0` (H5, and
`v_a ≠ v_b` because `(a,b)` is an `M′`-edge), `det(d,dS) = 0 ⟺ dS ∥ d`, so the three cases are
exhaustive and mutually exclusive, and `dS ∥ d, dS ≠ 0 ⟹ ⟨d,dS⟩ ≠ 0 ⟹ A ≠ 0`, which is exactly
the side-condition `A ≠ 0` that check.md R3.1a demanded of (T5.1e′). So **class 3 as reached from a
coincident-copy pair automatically satisfies R3.1a's side-condition** — a point the statement does
not make but which is the reason (T5.1e′) is safe here. Class 2's second root
`τ* = −B/A` follows from `g(τ) = τ(B + Aτ)`. **PROVED.**

**Cross-check against T5.3 (independent confirmation).** For a split pair,
`dS = S_{(v,g)} − S_{(v,f)} = J(2u_g − σ_g x_v) − J(2u_f − σ_f x_v) = 2 J Δu` (`σ_f = σ_g`).
Then by (K1) `A = σ_f⟨d, 2JΔu⟩ = −2σ_f det(d, Δu) = 2p_{T5.3}` ✓ and
`B = det(d, 2JΔu) = 2⟨d, Δu⟩ = 2r_{T5.3}` ✓ against (T5.2). The two independent derivations agree,
including both signs. **PROVED.**

**(L1.2b).** Type (b): `σ_g = −σ_f`, so `dS = J(2Δu − (σ_g−σ_f)x_v) = J(2σ_g x_{src(e)} − 2σ_g x_v)
= 2σ_g J(x_{src(e)} − x_v)` using (T1.3) `Δu = σ_g x_{src(e)}`; and by (K1)
`det(d, dS) = 2σ_g ⟨x_{src(e)} − x_{dst(e)}, d⟩`, so class 3 `⟺` the tested edge ⟂ the hinge edge.
**PROVED.** Type (c): `⟨d, Δu⟩ = 0`, which is T5.3's `r = 0` and core.md T5.H item 4 ("the flat
root is double"). **PROVED.** Type (d): the general formula
`dS = J(2Δu − (σ_g − σ_f)x_v)` is what I derived above before specializing. **PROVED.**

**(L1.2a) — I attacked this as the spec asked, and the attack fails.** Class 3 gives
`q = −p`, `r = 0`, hence `h(θ) = p(1 − cos θ)`; `h(0)=p+q=0`, `h′(0)=r=0`, `h″(0) = −q = p`, so the
statement's `p = h″(0)` is right. `1 − cos θ > 0` on `(0, 2π)`, so `sign h ≡ sign p ≠ 0` there,
including at `θ = π` (`h(π) = 2p`), which matters because H6 runs the deployment on `(0, π]`.
*The "never a contact" reading is correct and is not weakened by a permanent touching:* a contact
of the pair `(w,(a,b))` requires `w` to be collinear with `(a,b)`, i.e. `h_π(θ) = 0`; a *permanent*
touching would require `h_π ≡ 0`, which is the `A = B = C = 0` class, **not** class 3 (class 3 has
`A ≠ 0`). The two interval predicates of `contact.hpp` are only ever evaluated *at* a root, so no
root ⟹ no candidate angle ⟹ the deflated atom `true` is correct. **PROVED.** This also reproduces
Lemma T5.1e and its R3.1a side-condition, independently.

**(L1.2c).** `det(d, dS)` is quadratic in `t` (both factors affine), so `{det(d,dS) = 0}` is a
proper algebraic subset *unless the polynomial is identically zero* — and the whole point of the
sentence is that on the reference tilings it **is** identically zero (orthogonality forced by
symmetry). The statement says exactly this, so it is self-consistent, but "codimension ≥ 1" is only
justified for the pairs where the polynomial is not identically zero, and the statement does not
separate those two populations the way L1.1(iii)/(iv) does for `h(0)`. The count `1 138 of
2 145 387` is core.md's measurement, not re-run by me. **PLAUSIBLE-UNPROVED** (the qualitative
claim), **NOT-CHECKED** (the count).

## L1.3 — verdict 1⇔2⇔3 PROVED-INDEPENDENTLY · the boxed "so … are exactly" summary is **OVERSTATED (my id P1-L1-d)** · the split branch is **incomplete (P1-L1-e)**

**2 ⟹ 1.** `det(λ(Y_b−Y_a), Y_b−Y_a) = 0`. Trivial.

**1 ⟹ 2 (the only non-trivial implication, and it needs an argument the statement does not
supply).** `Y_b − Y_a = c C_ab + s S_ab = c d − σ_f s J d = R_f d ≠ 0` for every `θ`, since `R_f` is
a rotation and `d ≠ 0` (H5). So `h ≡ 0` gives `Y_w − Y_a = λ(θ)(Y_b − Y_a)` with an *a priori*
`θ`-dependent `λ`. Evaluate at the two ends: at `θ = 0`, `(c,s) = (1,0)` gives `C_aw = λ(0) d`; at
`θ = π`, `(c,s) = (0,1)` gives `S_aw = −σ_f λ(π) J d`. Substituting both back into
`Y_w − Y_a = c C_aw + s S_aw` and expanding in the basis `(d, Jd)` of `R²`, the vanishing of
`det(Y_b − Y_a, Y_w − Y_a)` at a general `θ` reads

```
        det( c d − σ_f s Jd ,  c λ(0) d − σ_f s λ(π) Jd )  =  −σ_f cs ( λ(π) − λ(0) ) det(d, Jd) = 0 .
```

`det(d, Jd) = |d|² ≠ 0`, so for any `θ` with `cs ≠ 0` we get `λ(π) = λ(0) =: λ`, and then
`Y_w − Y_a = λ(Y_b − Y_a)` identically in `θ`. **PROVED.**

**2 ⟺ 3.** With `λ` constant, 2 is `C_aw = λ d` **and** `S_aw = λ S_ab`. The first is exactly
3's collinearity. For the second, `S_aw = J(2Δu − σ_g x_{v_w} + σ_f x_{v_a})` and
`λ S_ab = −λ σ_f J d = −σ_f J(x_{v_w} − x_{v_a})` (using the first). `J` is invertible, so
`2Δu − σ_g x_{v_w} + σ_f x_{v_a} = −σ_f x_{v_w} + σ_f x_{v_a}`, i.e.
`2Δu = (σ_g − σ_f) x_{v_w}`. **Exactly 3. PROVED.**

**The hinge specialization (L1.9).** `σ_g = −σ_f` and `Δu = σ_g x_{src(e)}` turn 3's second
condition into `2σ_g x_{src(e)} = 2σ_g x_{v_w}`, i.e. `x_{v_w} = x_{src(e)}` (as points; as
*vertices* this uses H1, the embedding being injective). Together with 3's collinearity this is
(L1.9). **PROVED.**

**The split specialization — incomplete (P1-L1-e).** `σ_g = σ_f` turns 3's second condition into
`Δu = 0`, correctly; but the statement prints "it reads `Δu = 0`", dropping 3's **first** condition.
`Δu = 0` alone does **not** give `h ≡ 0`: one still needs `x_{v_w}` on the line of `(x_{v_a},
x_{v_b})`. Compare the hinge branch, which does print both halves. This is a **gap in the printed
statement**, not an error in the equivalence.

**The boxed summary is overstated (P1-L1-d).** The box asserts the identically-zero pairs are
*exactly* the hinge-point pairs. But condition 3 is satisfied by at least two further families that
the box excludes:
* `f, g` split-adjacent (or any two faces with `σ_f = σ_g`) with `Δu = 0` **and** `v_w` collinear
  with `(v_a,v_b)`. `Δu = 0` means the split cut never opens — `dS_e = 2JΔu = 0`, exactly the
  `q_e = 0` degenerate case `zero_plus.hpp` is written to detect — so it is a real, describable
  configuration, not a vacuous one.
* type (d), `σ_g = −σ_f`, two faces meeting only at `v_w`, with `2Δu = 2σ_g x_{v_w}` and the
  collinearity. Nothing in H1–H6 forbids it.
Neither family is generic, but "exactly" is a claim of *equality of sets*, and the equivalence
1⇔2⇔3 does not deliver it. Either the box needs a genericity hypothesis (`Δu ≠ 0` on every
split-adjacent pair, i.e. every split cut opens — which is precisely what K5/F30 says **fails** in
practice), or it must be widened to include the `Δu = 0` families. I grade this
**PLAUSIBLE-UNPROVED as printed / PROVED once widened**; the equivalence itself is PROVED. My test
T-L1-b below searches the corpus for a witness.

**The spec's suggested attack** — "find a pair with `p = q = r ≡ 0` that is not a listed permanent
incidence" — is therefore answered *structurally* by P1-L1-d: the candidates are the `Δu = 0` split
pairs and the type-(d) pairs, and whether they occur is a corpus question (T-L1-b).

## L1.4 — verdict PROVED-INDEPENDENTLY, conditional on L1.1(ii) and L1.3

(i) is a restatement: the class is decided by `C = h(0)` and `B = 2h′(0)` (both read off `(p,q,r)`),
and on coincident-copy pairs (L1.8) makes both functions of `dS` alone. The three bullets follow
from L1.3, L1.1(ii) and L1.2 respectively, so (i) inherits P1-L1-a/d. (ii) is then immediate. (iii)
is an honest caveat and I agree with its *direction* argument: routing a `C ≡ 0` pair to (T5.1c)
makes `g(0)·g(T) > 0` fail, so the pair is reported as having a root, so `NOROOT` is refused, so
`R(ε)` shrinks — conservative, i.e. soundness-preserving. This matches core.md T5.2b.2 caveat (1)
and check.md R3.1a. **PROVED** (modulo the two gaps above).

## L2.1 — verdict PROVED-INDEPENDENTLY, with one hypothesis that must be stated (P1-L2-a)

**Derivation of (L2.8), independently.** Under a shape-space perturbation, `P₀` is *fixed data*
(the quotient system fixes the lattice `T`), so `δK = δQ · P₀⁻¹`. By (L2.4) `w_τ` is affine in `X`
with linear part `Σ_{path} σ_{head(e)} x_{src(e)}`, applied to each coordinate separately. The
design coordinate `t_{2i}` moves `X` by `φ_i e_xᵀ`, so `δw_τ = d_τ(φ_i) e_x`; `t_{2i+1}` moves `X`
by `φ_i e_yᵀ`, so `δw_τ = d_τ(φ_i) e_y`. With `Q = 2J[w_h w_v]`,
```
   ∂Q/∂t_{2i}   = 2 J e_x d_iᵀ = 2 e_y d_iᵀ           (J e_x = e_y)
   ∂Q/∂t_{2i+1} = 2 J e_y d_iᵀ = −2 e_x d_iᵀ          (J e_y = −e_x)
```
and right-multiplying by `P₀⁻¹` gives **exactly (L2.8)**, both formulas and both signs. **PROVED.**
This also matches the code: `kill_b3.cpp::achievable` reads `D(i,0..1) = ½ (M_{2i} T)(row 2)`,
which is `½ · 2 d_iᵀ = d_iᵀ` — the same convention.

**(L2.9), and the spec's attack on the factorisation.** Put `g_iᵀ := d_iᵀ P₀⁻¹`. Then
`M_{2i} = 2 e_y g_iᵀ`, `M_{2i+1} = −2 e_x g_iᵀ`, so
```
   span{M_j} = { e_x aᵀ + e_y bᵀ : a, b ∈ G } ,      G := span{g_i} ⊆ R² ,
```
because the `x`-row and the `y`-row of the generated matrices can be chosen independently: the
`t_{2i+1}` generators populate the first row only and the `t_{2i}` generators the second row only.
Hence `rank(A) = dim span{M_j} = 2 dim G`, and `dim G = rank(D P₀⁻¹) = rank(D)` by invertibility of
`P₀` (H-L2). So `rank(A) = 2 rank(D)`, and `dim 𝒦 ∈ {0,2,4}`. **PROVED.**
*The attack fails.* The `x`-copy and the `y`-copy of one null vector give `−2 e_x g_iᵀ` and
`2 e_y g_iᵀ`, which for `g_i ≠ 0` are linearly independent **for every `P₀`** — they live in
complementary rows. A square lattice, or any other special `P₀`, only changes `g_i = P₀⁻ᵀ d_i`,
never the `e_x`/`e_y` split, so the two cannot collapse. Independence fails only when `g_i = 0`,
i.e. `d_i = 0`, in which case *both* vanish and the pair contributes `0`, not `1` — which is why
the rank is always **even**, and why `dim 𝒦 = 1` or `3` is impossible.

**P1-L2-a — a hypothesis the statement needs and does not print.** (L2.5)'s path-independence
argument is stated as "a change of path adds a cycle of `Γ`, and the closure of `δ` on cycles is
exactly `L φ = 0`". On the **torus quotient** this is false as written: `Γ` has non-contractible
cycles, whose `δ`-sums are precisely the non-zero `w_τ` the lemma is about, and those are **not**
rows of `L` (`quotient_system` builds `L` from the `H` hole rows, with the lattice offsets carried
on the right-hand side). The correct statement is that two walks `f₀ → f₀+τ` differ by a closed walk
of **zero total lattice offset**, and *those* are generated by the hole cycles, hence killed by
`L φ = 0`. With that amendment the argument is right; as printed it proves too much. **gap.**

## L2.2 — verdict (a) PROVED (trivial, and honestly so) · (b) PROVED, but only via a fact the statement omits (P1-L2-b) · (c) PROVED, same caveat · (d) NOT-CHECKED (measurement; my test T-L2-d re-measures it) · (e) STATEMENT-ILL-POSED as a claim, acceptable as a labelled conjecture

**(a).** `D ∈ R^{k×2}` so `rank(D) ≤ min(2,k)` **by shape** — the spec is right that this is
trivial. What (a) actually adds is the *transport* through L2.1: `dim 𝒦 = 2 rank(D) ≤ 4` and
`≤ 2k`. Both are correct consequences. "Both bounds are attained" is a measurement. **PROVED**
(as a consequence), **NOT-CHECKED** (attainment).

**(b) and P1-L2-b.** `rank(D) = 0 ⟺` every entry `d_τ(φ_i) = 0 ⟺` (by linearity) `d_h` and `d_v`
vanish on `span{φ_i}`. But `span{φ_i}` is `ker[L; e_pin]`, **not** `ker L`, and the statement says
"for every `φ ∈ ker L`". These differ by exactly one dimension: constants lie in `ker L` (every hole
row of `L` has coefficients summing to zero, Eq. (2) being a sum of edge differences) and the pin
row removes them. So (b) as printed is a *strictly stronger* right-hand side than the left-hand
side supports — **unless** `d_τ(1) = 0`. It is: along any walk `f₀ → f₀+τ` in `Γ`, `σ` alternates
(§0.5, `Γ` bipartite), and `f₀+τ` is the same face class as `f₀`, hence carries the same `σ`, hence
every such walk has **even length**, hence `Σ_{path} σ_{head(e)} = 0` and `d_τ(1) = 0`. So the
equivalence survives, but it rests on a parity argument the statement does not print. **PROVED,
with the missing step supplied.** (Geometrically this is just "a global translation does not change
`K`", which it had better not.)

**(c).** Let `R = row(L) ⊆ (R^{n_q})*`. Since `(row L)^⊥ = ker L`, a covector vanishes on `ker L`
iff it lies in `R`. Hence
`rank[L; d_h; d_v] − rank(L) = dim((R + span{d_h,d_v})/R) = dim span{ d_h|_{ker L}, d_v|_{ker L} }`,
which is the rank of the map `ker L → R²`, `φ ↦ (d_h(φ), d_v(φ))`. And `rank(D)` is the rank of the
**same** map restricted to `ker[L;e_pin]`. The two agree because `ker L = ker[L;e_pin] ⊕ span{1}`
and `d_τ(1) = 0` by the parity argument of P1-L2-b. **PROVED**, again modulo that step — without it
(c) could be off by one. The reading "`rank(D)` counts how many period covectors are not already
determined by Eq. (2)" is then exactly right.

**(d).** A measurement plus an interpretation. The interpretation is sound: `rank(D) = 0` gives
`dim 𝒦 = 0` by L2.1 regardless of `k`, so `min(4, 2·dim_null)` is refuted by any `rank(D) = 0`
pattern with `k > 0`. `tr K₀ = 0, det K₀ = 1` forces `K₀` to satisfy `K₀² = −I` (Cayley–Hamilton),
and a real `2×2` with `K₀² = −I` and `det = +1` is conjugate to `J`; combined with
`J(θ) = cI + sK₀` this is a rotation by `θ/2` **iff `K₀ = J` exactly**, not merely `K₀² = −I`
(e.g. `K₀ = −J` also has `tr = 0`, `det = 1`, and gives rotation by `−θ/2`). The printed two
numbers do not distinguish `J` from `−J`; the conclusion "rigid rotation by `θ/2`" therefore needs
the sign of `K₀₂₁`, which is not printed. **Minor gap, P1-L2-c**; the load-bearing claim
(`dim 𝒦 = 0` with `k = 4`) is unaffected. My test T-L2-d checks `rank(D)`, `K₀` and the sign.

**(e).** As a mathematical statement, "`rank(D) = min(2,k)` for every pattern **except those whose
symmetry forces `D = 0`**" is **ILL-POSED**: "symmetry forces `D = 0`" is not a predicate on
`(M, σ, T)` — read literally the exception clause is "except those where it fails", which makes the
statement a tautology. It is correctly *labelled* a conjecture and correctly reports that the author
could not find the separating predicate, and the evidence (132/133) is a fact worth printing. So:
**STATEMENT-ILL-POSED as a claim; acceptable as printed provided it is never used.** My counting
test T-L2-e re-measures the 132/133 split.

## L2.3 — verdict PROVED-INDEPENDENTLY, with one degenerate case unstated

The similarities of `R²` are `span{I, J}` (as a linear space; the *invertible* ones are that space
minus `0`). `J(θ) = cI + sK`. If `K ∈ span{I,J}` then so is `J(θ)` for every `θ`. Conversely if
`J(θ₀) ∈ span{I,J}` for one `θ₀` with `s ≠ 0`, then `K = (J(θ₀) − cI)/s ∈ span{I,J}`. And
`K ∈ span{I,J} ⟺ K₁₁ = K₂₂ ∧ K₁₂ = −K₂₁`, two **linear** equations in `t` because `K` is affine in
`t` (L2.1). **PROVED.** Unstated: a *conformal* map is usually required to be invertible, and
`cI + sK` can be **singular** at an interior `θ` even with `K` a similarity — e.g. `K = −I` gives
`J(π/2) = 0`. So "similarity for every `θ`" is right for the linear-span reading and needs
`det J(θ) ≠ 0` added for the conformal reading. **(wording, P1-L2-d.)**

## Pass-1 summary table

| statement | blind verdict |
|---|---|
| L1.1(i) | PROVED-INDEPENDENTLY |
| L1.1(ii) implication | PROVED-INDEPENDENTLY |
| L1.1(ii) "exactly four types" | PLAUSIBLE-UNPROVED — gap P1-L1-a (type (d) must say "no edge *containing v_w*") |
| L1.1(iii) | PROVED-INDEPENDENTLY (H-L1 not needed for the stated dichotomy; P1-L1-b) |
| L1.1(iv) | NOT-CHECKED (measurement); arithmetic internally consistent; P1-L1-c |
| L1.2 (L1.8) + trichotomy | PROVED-INDEPENDENTLY (both endpoint cases; cross-checked against T5.3) |
| L1.2a | PROVED-INDEPENDENTLY; the "never a contact" attack fails |
| L1.2b | PROVED-INDEPENDENTLY (all three types) |
| L1.2c | PLAUSIBLE-UNPROVED (qualitative); count NOT-CHECKED |
| L1.3 (1⇔2⇔3) | PROVED-INDEPENDENTLY |
| L1.3 split branch | gap P1-L1-e (collinearity dropped) |
| L1.3 boxed "exactly" | PLAUSIBLE-UNPROVED — P1-L1-d (`Δu = 0` split pairs and type (d) pairs also satisfy 3) |
| L1.4 | PROVED-INDEPENDENTLY, inheriting P1-L1-a/d |
| L2.1 (L2.8)+(L2.9) | PROVED-INDEPENDENTLY; the independence attack fails for every `P₀` |
| L2.0c path-independence | gap P1-L2-a (needs "zero lattice offset", not "any cycle of `Γ`") |
| L2.2(a) | PROVED (trivial by shape, as the spec suspected) |
| L2.2(b) | PROVED via the omitted parity step P1-L2-b (`ker L` vs `ker[L;e_pin]`) |
| L2.2(c) | PROVED, same caveat |
| L2.2(d) | NOT-CHECKED; minor gap P1-L2-c (`K₀ = J` vs `−J` not pinned by `tr`/`det`) |
| L2.2(e) | STATEMENT-ILL-POSED as a claim; acceptable as a labelled conjecture |
| L2.3 | PROVED-INDEPENDENTLY; P1-L2-d (invertibility) |

**No COUNTEREXAMPLE was found to any statement in Pass 1.**

---

# PASS 2 — comparison against `derivations/lemmas.md`

Opened only after the Pass-1 section above was complete on disk.

## Statements-file hygiene

`derivations/lemmas_statements.md` is a **verbatim** subset of `derivations/lemmas.md`: I diffed
every non-blank line of every extracted block against the full file mechanically, and the only line
that is not present verbatim is the orchestrator's own extraction note in the header. **No statement
discrepancy.** One cosmetic point: the extracted L1.0 heading reads "one orthogonal frame kills all
four cases", but the frame itself ([A] L1.0b–d) is withheld, so the heading advertises material the
blind reader cannot see. That is what made the blind pass genuinely blind, and it is why my Pass-1
proof of L1.3 `(1) ⟹ (2)` is a clumsier two-endpoint argument than the Deriver's.

## L1.0 — the orthogonal frame (withheld from Pass 1)

| step | verdict | note |
|---|---|---|
| L1.0b `C_ab = d`, `S_ab = −σ_f J d`, `D₀ = −σ_f\|d\|² ≠ 0` | **AGREE** | T3.3 + H5; the identity needs `a,b` in one face, which the Deriver flags himself as weak step 2 |
| L1.0c `C = βD₀`, `A = −γD₀`, `B = (δ−α)D₀` | **AGREE** | re-derived: `det(C_ab,αC_ab+βS_ab) = βD₀`; `det(S_ab,γC_ab+δS_ab) = −γD₀`; `det(C_ab,S_aw)+det(S_ab,C_aw) = δD₀ − αD₀` |
| L1.5 the four-way linear split | **AGREE** | `D₀ ≠ 0`, so the class conditions are linear in `(α,β,γ,δ)` |
| L1.0d `h′(0) = r`, `h″(0) = −q = (A−C)/2` | **AGREE** | and the class split is chart-free, which does remove check.md D10 from the classification |

This frame is a genuine improvement on my Pass-1 route and I withdraw nothing from Pass 1 because of
it — every Pass-1 conclusion is reproduced by it.

## L1.1 — step-by-step

| proof step | verdict | note |
|---|---|---|
| 1 (`h(0) = det` of source vertices) | **AGREE** | identical to my Pass-1 step |
| 2 (coincident copy ⟹ `h(0) = 0`) | **AGREE** | no hypothesis on `X`, as claimed |
| 3 (the four types are exhaustive) | **AGREE with the proof, DISAGREE with the printed table** — **D-L1-1** | see below |
| 4 (degree ≤ 2, measure zero) | **AGREE** | |
| 5 ((iv) is `[N]`) | not re-run by the Deriver's program; re-measured independently, see **D-L1-4** | |

**D-L1-1 — wording, L1.1(ii) type (d).** Proof step 3 states the correct dichotomy: two faces
meeting at `v_w` "either share an edge **at that vertex** … or meet only at the vertex". The
statement's table row (d) instead reads "`f, g` share the vertex `v_w` but **no edge**", which does
not cover `f, g` adjacent across an edge `e` while also meeting at some `v_w ∉ e`. This is my
Pass-1 finding P1-L1-a and it is a *wording* defect only: the proof is right, one word fixes the
table. **Measured 0 occurrences** (test `L1-k3`, 835 272 coincident-copy pairs over 60 graphs).

**D-L1-6 — wording, L1.1(iii).** Headed "(converse, under H-L1)", but the dichotomy stated is
unconditional; H-L1 is the separate assertion that the first branch is empty. Pass-1 P1-L1-b.

## L1.2 — step-by-step

| proof step | verdict | note |
|---|---|---|
| 1 (`dC = 0`, `Y_w − Y_{a*} = s·dS`) | **AGREE** | matches `zero_plus.hpp`'s `dC_e = 0`, correctly generalised |
| 2 (`α,β` in both endpoint cases; `A = −γD₀`, `B = δ′D₀`) | **AGREE** | I checked the `v = v_b` case separately: `δ = 1 + δ′`, `α = 1`, so `B = (δ−α)D₀ = δ′D₀`. This is the one place a sign could hide and it does not |
| 3 (`r = ½det(d,dS)`) | **AGREE** | |
| 4 (`p = ½σ_f⟨d,dS⟩` via `\|d\|² = −σ_f D₀`) | **AGREE** | identical to my Pass-1 result, reached differently |
| 5 (trichotomy; `A ≠ 0` free) | **AGREE**, and this is the best step in Part L1 | it turns check.md R3.1a's hand-added side-condition into a consequence: `A = B = 0 ⟺ dS = 0` |
| 6 (`τ* = −B/A`, and T5.4) | **AGREE** | I verified the T5.3 cross-check independently in Pass 1: `A = 2p_{T5.3}`, `B = 2r_{T5.3}` |
| 7 ((L1.2a)) | **AGREE** | reproduces Lemma T5.1e; my Pass-1 attack on the "never a contact" reading failed |
| 8 ((L1.2b) types (b) and (c)) | **AGREE** | both closed forms re-derived independently in Pass 1 |
| 9 (`[N]`) | re-measured independently, see Tests | |

**D-L1-5 — I agree with the Deriver against `core.md`.** L1.2's "What it replaces" calls
`core.md` T5.2b.2's sentence — "class 3 is the sub-case in which the flat configuration additionally
satisfies `r ≡ 0`, which by (T3.2) is again a statement about which `M′`-copies coincide, **not about
`t`**" — wrong as printed. It is wrong: `r = ½det(d,dS)` is a quadratic polynomial in `t`. My run
independently finds class-3 pairs that appear at some `t` and not others, and finds 1 185 class-3
occurrences that are not coincident-copy pairs at all. Not a disagreement with `lemmas.md`; recorded
because it is a correction `core.md` must absorb.

## L1.3 — step-by-step

| proof step | verdict | note |
|---|---|---|
| 1 (`(1) ⟹ (2)` via the frame) | **AGREE**, and it is cleaner than my Pass-1 argument | `β = γ = 0, δ = α` gives `C_aw = αC_ab`, `S_aw = αS_ab` with the *same* `α`, hence a constant `λ`. The parenthetical about orthogonality doing the work is correct |
| 2 (`(2) ⟹ (1)`) | **AGREE** | |
| 3 (`(2) ⟺ (3)`) | **AGREE** | identical to my Pass-1 computation, including the cancellation of `σ_f x_{v_a}` |
| 4 (hinge specialization, (L1.9)) | **AGREE** | |
| 5 (split specialization) | **AGREE with the algebra, DISAGREE with the conclusion drawn** — **D-L1-3** | |
| 6 (non-adjacent) | **AGREE with the algebra, DISAGREE with the conclusion drawn** — **D-L1-3** | |
| 7 (`[N]`) | re-measured independently; both directions reproduced exactly | |

**D-L1-3 — gap, and an internal contradiction. The strongest disagreement in Part L1.**
The boxed summary asserts the identically-zero pairs are **exactly** the hinge-`src` pairs. The
Deriver's own proof steps 5 and 6 concede two further families that satisfy condition 3:

* step 5: `f, g` split-adjacent with `Δu = 0`. The step does not argue that this cannot happen — it
  argues that when it happens `Θ_max = 0` and "such a design is rejected by the very first predicate
  of the certificate". That is a statement about the *certificate*, not about whether the pair has
  `h ≡ 0`. It therefore does not support "exactly".
* step 6: `f, g` non-adjacent, where `2Δu = (σ_g−σ_f)x_{v_w}` "can hold accidentally at isolated
  `t`". Again a concession, not an exclusion.

So the box contradicts steps 5–6. What steps 5–7 actually prove is: *on this corpus*, measured
`0` split-adjacent and `0` non-adjacent instances. Two further consequences:

1. The box's closing sentence — "This is a **combinatorial** description: decided by `(M, σ)` and
   the flat incidences, **not by `t`**, except for the single collinearity in (L1.9)" — is false for
   the split branch, because `Δu = u_g − u_f` is *linear in `X`*, so `Δu = 0` is a second condition
   on `t`, not a combinatorial one.
2. The clean repair is a named hypothesis, and the natural one is **"every split cut opens",
   `Δu ≠ 0` on every split-adjacent face pair** — which is exactly the condition
   `zero_plus.hpp`/K5/F30 report as *failing* on measured designs. So the hypothesis is not free and
   should be printed rather than absorbed.

My tests confirm the corpus claim in both directions with zero exceptions (`L1-i`, `L1-j`:
223 848 / 223 848), so the practical content of L1.3 stands. What does not stand is the word
"exactly" as a theorem.

**D-L1-2 — wording, L1.3's split branch.** The statement prints "if `f` and `g` are
**split-adjacent**, it reads `Δu = 0`", dropping condition 3's collinearity half, which the hinge
branch does print in full. Proof step 5 says "condition 3's **second half**", so the proof is
correct and only the statement is elliptical. Pass-1 P1-L1-e.

## L1.4 — **AGREE** throughout, inheriting D-L1-1 and D-L1-3

(i) is L1.1(i)+(ii), L1.2 step 5 and L1.3, as claimed; (ii) follows; (iii)'s conservativeness
argument is correct — `C ≡ 0` makes (T5.1c)'s first atom `g(0)·g(T) > 0` false, so `NOROOT` is
refused and `R(ε)` shrinks. The Deriver's own "weakest step 1" caveat about the *near-threshold*
population is fair and I did not close it either: my tests classify with the same
`1e−11·scale²` rule and therefore inherit it.

## L2.0 — step-by-step

| step | verdict | note |
|---|---|---|
| L2.0a (the period law, `u_{f+τ} = u_f + w_τ + σ_f τ/2`) | **AGREE** | I re-checked the ansatz on a hinge edge: required `σ_g(x_{src}+τ)`, ansatz gives `σ_g x_{src} + (σ_g−σ_f)τ/2 = σ_g x_{src} + σ_g τ`. Uniqueness up to a constant needs `Γ` connected, correctly cited |
| L2.0b (`P_θ = cP₀ + sQ`, `K = QP₀⁻¹`) | **AGREE** | the `σ_f τ` terms do cancel |
| L2.0c (L2.4)–(L2.6) | **AGREE with the formulas, DISAGREE with the path-independence justification** — **D-L2-1** | |
| L2.0d (L2.7), coordinate separation | **AGREE** | and this is the Deriver's own "weakest step 3"; I closed it, see **D-L2-5** |

**D-L2-1 — gap, L2.0c.** "its restriction to `ker L` is independent of the path … because a change
of path adds a cycle of `Γ`, and the closure of `δ` on cycles is exactly `L φ = 0`". On the **torus
quotient** this proves too much: `Γ` carries cycles of non-zero lattice offset whose `δ`-sum is
precisely the non-zero `w_τ` the lemma is about, and those are **not** rows of `L` — `quotient_system`
puts the offsets on the right-hand side. The correct statement is that two walks `f₀ → f₀+τ` differ
by a closed walk of **zero total lattice offset**, and those are generated by the hole cycles.
Pass-1 P1-L2-a, and now measured: on **120 of 134** patterns the two path covectors differ as
covectors on `R^{n_q}`, while their restrictions to the shape space agree to `1.37e−15` (test
`L2-h`). The conclusion holds; the printed reason does not.

## L2.1 — step-by-step: **AGREE on every step**

| proof step | verdict |
|---|---|
| 1 (`K` affine; `P₀` fixed) | **AGREE** |
| 2 (`M_{2i}P₀ = 2Je_x d_iᵀ = 2e_y d_iᵀ`) | **AGREE** — identical to my Pass-1 derivation |
| 3 (`M_{2i+1}P₀ = 2Je_y d_iᵀ = −2e_x d_iᵀ`) | **AGREE** |
| 4 (the span, `(a,b)` ranges over `row(D) × row(D)`) | **AGREE** |
| 5 (injectivity of `Ψ`; where `P₀` singular breaks it) | **AGREE** — and this is exactly the attack the spec asked me to run; it fails, because the two generators occupy **complementary rows** of the matrix, so no `P₀` can collapse them |
| 6 (`dim 𝒦 = 2 rank(D)`, parity) | **AGREE** |
| 7 (`[N]`) | see **D-L2-5** |

Lemma L2.1 is the strongest result in the file: correct, complete, and it yields more than the
dimension formula it replaces.

## L2.2 — step-by-step

| proof step | verdict | note |
|---|---|---|
| 1 ((a)) | **AGREE** | trivial by shape, and the file says so |
| 2 ((b)) | **AGREE on the conclusion; the inference as printed is WRONG** — **D-L2-2** | |
| 3 ((c)) | **AGREE on the conclusion; one printed step is FALSE** — **D-L2-2** | |
| 4 ((d) `[N]`) | re-measured; **D-L2-3** | |
| (e) | **STATEMENT-ILL-POSED** as a claim — **D-L2-4**; correctly labelled a conjecture | |

**D-L2-2 — a false step, in both (b) and (c). The strongest disagreement in Part L2 on the
mathematics.** Proof step 3 writes "Let `N := ker L`… the `φ_i` are a basis of `N`". They are not.
`Φ` is, by L2.0's own set-up sentence, "a null-space basis of the **scalar** system `[L; e_pin]`",
and `ker[L;e_pin] ⊊ ker L`: the constant vector `1` lies in `ker L` (every hole row of `L` has
coefficients summing to zero, Eq. (2) being a sum of edge differences) and is removed by the pin. So
`dim ker L = dim span{φ_i} + 1` and step 3's identification is off by one dimension. The same
conflation is in step 2 of (b): "for every basis vector, hence by linearity for every `φ ∈ ker L`".

Both conclusions **survive**, but only because of a fact `lemmas.md` never states:
```
        d_τ(1) = 0 .
```
*Proof.* `Γ` is bipartite and `σ` is a proper 2-colouring (core.md 0.5). `f₀+τ` is the same face
class as `f₀`, hence carries the same `σ`, hence every lifted walk `f₀ → f₀+τ` has **even length**,
hence its alternating sum `Σ σ_{head(e)}` is zero. ∎ Geometrically: a global translation of the
pattern must not change `K`, and this is why it does not. Pass-1 P1-L2-b, now **measured exactly
`0.000e+00` on 134 patterns** (test `L2-g`). Two lines in `lemmas.md` close it.

**D-L2-3 — gap, minor, L2.2(d).** "`tr K₀ = 0.000000` and `det K₀ = 1.000000`, i.e. `K₀` is a
rotation by `π/2`". Those two numbers give `K₀² = −I` (Cayley–Hamilton), which is satisfied by
`+J` **and** by `−J`; the latter would make `J(θ)` a rotation by `−θ/2`. The sign of `K₀(2,1)` is
what decides, and it is not printed. Measured: `K₀(1,0) = +1.000000`, so the claim is **true** — it
just does not follow from the printed evidence. Pass-1 P1-L2-c.
*Second remark on the same worked case:* `squares_3x3`'s own `K₀` is measured on a `3×3` super patch
whose hinge graph has **12 components**, with the `(0,0)` and `(1,0)` copies of a face in different
ones. Its `PeriodicJac::consistency` is nonetheless exactly `0`, so the measurement is sound there —
but the pattern L2.2(d) is built on is one super-patch step away from the defect of **D-L2-5**.

**D-L2-4 — ill-posed, L2.2(e).** "for every pattern **except those whose symmetry forces `D = 0`**"
is not a predicate on `(M, σ, T)`; read literally the exception clause says "except where it fails",
which makes the sentence a tautology rather than a conjecture. The *evidence* is real and I
reproduce it (142/143, sole exception `squares_3x3`). It is correctly labelled CONJECTURE and
correctly not used anywhere, so this is a defect of phrasing, not of honesty. Pass-1 P1-L2-e.

**D-L2-5 — substantive, and it closes the Deriver's own weakest step 3.**
`lemmas.md`'s "For the Checker / Weakest step 3" says, of recomputing `D` directly from (L2.5):
"I did **not** do that; my `D` is the one the code computes." I did it (test `L2-f`; the routine is
`dtau_covector`, which walks the **lifted** hinge graph from `(f₀, 0)` to `(f₀, τ)` and sums
`σ_{head}` at the source class, never touching `K`, `M_j` or `AchievableSet`).

* On **132 of 134** patterns with `k > 0` the independently rebuilt `D` matches
  `AchievableSet::D` to **`2.66e−13`** relative. (L2.7), (L2.8) and the reading of `D` off the
  **even** generators are therefore confirmed **non-circularly**. This is the single most useful
  thing in my report for Part L2: the lemma's factor of two is now attacked and survives.
* On **2 patterns the two disagree by `3.72e−2` and `1.63e−2`** relative. One of them,
  **`snub_square_3x3`, is a K7 pattern — i.e. it is inside Check L2's own "133 / 133"**.

The cause is not the lemma. `PeriodicJac::consistency` — the spread of the measured period over
`(face, corner)`, and the header's own "numerical witness that the deployed tiling really is a
translate of itself" — is `2.750e+01` and `1.078e+01` on those two, against `≤ 6.4e−14` on every
other pattern in the corpus. The mechanism, traced: `build_super(q, 1)` produces a `3×3` patch whose
hinge graph has **9** (resp. **4**) components, and the `(0,0)` and `(1,0)` copies of a cell face
land in **different** components, so their face potentials carry **independent additive constants**
— exactly core.md **T1.H.1**. `H2` (`Γ` connected), a *standing hypothesis of L2.1*, therefore fails
on the object `K` is measured on, even though the **quotient** `Γ_q` is connected (I checked: 1
component, quotient residual `2.7e−15`).

`achievable()` in `kill_k7.cpp` / `kill_b3.cpp` **never reads `consistency`**. And note why the
`rank(A) = 2 rank(D)` equality does not catch it: both sides are read off the same `M_j`, so the
equality is insensitive to whether those `M_j` mean anything. The `(L2.8)` row-vanishing test
(`1.59e−14`) is likewise insensitive. Only a `D` built from a different definition sees it.

Concrete consequences, in the order they matter:
1. **Check L2's "`rank(A) = 2 rank(D)` on 133 / 133" should read 132 / 133 measured + 1 excluded**,
   and Lemma L2.1 should carry the sentence "`H2` must hold on the patch the period is measured on,
   not merely on the quotient; `PeriodicJac::consistency` is the test."
2. `achievable()` should be gated on `consistency ≤ tol`, or `build_super` given a larger `half`
   until `Γ_super`'s components no longer separate a face from its translates.
3. `squares_3x3`, the pattern L2.2(d)/(e) both turn on, is in the same family of super-patch
   disconnection (12 components) and survives only because its `D` is exactly `0`.
`c(Γ_super) > 1` on its own is common (**26 of 143**) and usually harmless; the sharp predicate is
`consistency`, and `n_Dind_bad == n_inconsistent` is asserted in the test.

## L2.3 — **AGREE on the statement and on both implications; one printed step is FALSE** — **D-L2-6**

The forward and converse implications are right, `similarities = span{I,J}` is right, and (L2.11)
is right and linear in `t` (`K` affine by L2.1 step 1).

**D-L2-6 — a false step.** The proof's non-singularity argument reads: "`det J(θ) = (c+sα)² +
(sβ)²`, which vanishes only if `c + sα = sβ = 0`; the map is then the zero matrix, excluded since
`det J(0) = 1` and `det J(θ)` is a non-negative quadratic in `(c,s)` **that vanishes at most at
isolated `θ`**." The final clause concedes exactly what the sentence sets out to exclude, and
`det J(0) = 1` does not prevent a vanishing at a later `θ`. Concretely, `K = −I` is a similarity
(`α = −1`, `β = 0`) and gives `J(π/2) = 0`, with `π/2` **inside** the deployment range `(0, π]` of
H6. So under the reading the proof itself adopts — similarities are `span{I,J}` *minus the origin* —
"`J(θ)` is a similarity for **every** `θ`" is false there. The repair is one of: state the
conclusion for the linear span (which is all L2.3 is used for), or add the side-condition
`K ∉ {αI : α < 0}`. I did **not** search the corpus for a pattern realising `K = αI` with `α < 0`,
so I do not claim the case occurs; the defect in the printed proof is independent of that.

## Pass-2 disagreement list

| id | kind | one line |
|---|---|---|
| **D-L1-1** | wording | L1.1(ii) row (d) should read "no edge **containing `v_w`**"; the proof already has it right. Measured 0 occurrences. |
| **D-L1-2** | wording | L1.3's split branch drops condition 3's collinearity half; the proof (step 5) has it. |
| **D-L1-3** | **gap** | L1.3's boxed "**exactly** the hinge-point pairs" is contradicted by its own proof steps 5–6, which concede split-adjacent `Δu = 0` and non-adjacent pairs; and "not by `t`" is false for the split branch, `Δu` being linear in `X`. Proved only *on the corpus* (which my tests confirm, 223 848/223 848 both ways). |
| **D-L1-4** | wording | L1.1(iv)'s counts do not state their unit. My independent run reproduces the accidental set **exactly** (1 116) and the rigid set closely (20 447 vs 22 251) as **pair** counts, but gets **7 624** structural pairs against the printed **174 354** — which reconciles only if that one number is an *instance* count. |
| **D-L1-5** | agreement | I confirm the Deriver's correction of `core.md` T5.2b.2 ("class 3 … not about `t`" is wrong; `r = ½det(d,dS)` is a polynomial in `t`). |
| **D-L1-6** | wording | L1.1(iii) is labelled "under H-L1" but states an unconditional dichotomy. |
| **D-L2-1** | **gap** | L2.0c's path-independence reason ("a change of path adds a cycle of `Γ`") is false on the torus; it needs "a closed walk of **zero lattice offset**". Measured: the two path covectors differ off the null space on 120/134 patterns and agree on it to `1.37e−15`. |
| **D-L2-2** | **wrong step** | L2.2(c) proof step 3 (and (b) step 2) treat `{φ_i}` as a basis of `ker L`; they are a basis of `ker[L;e_pin]`, one dimension smaller. Both conclusions survive **only** via the unstated `d_τ(1) = 0`, which follows from the bipartiteness of `Γ` and is measured exactly `0` on 134 patterns. |
| **D-L2-3** | gap, minor | L2.2(d): `tr K₀ = 0` and `det K₀ = 1` do not distinguish `K₀ = J` from `−J`. Measured `K₀(1,0) = +1`, so the claim is true but unsupported by the printed numbers. |
| **D-L2-4** | ill-posed | L2.2(e)'s exception clause "except those whose symmetry forces `D = 0`" is not a predicate. Evidence reproduced 142/143. |
| **D-L2-5** | **substantive** | The Deriver's weakest step 3, closed: `D` rebuilt from (L2.5) matches the code on **132/134** to `2.66e−13`, and **disagrees by `3.7e−2` on `snub_square_3x3`**, a pattern inside Check L2's own "133/133", because `H2` fails on the `3×3` super patch (`c(Γ_super) = 9`, `PeriodicJac::consistency = 27.5`) that `K` is measured on. `achievable()` never reads `consistency`, and `rank(A) = 2rank(D)` cannot detect this because both sides come from the same `M_j`. |
| **D-L2-6** | **wrong step** | L2.3's non-singularity argument concedes what it claims to exclude; `K = −I` gives `J(π/2) = 0`, inside H6's range. |

Nothing in Pass 1 is withdrawn in Pass 2. Every Pass-1 verdict is either confirmed by the Deriver's
proof (L1.1(i)(ii)(iii), L1.2 in full, L1.3's equivalence, L2.1 in full, L2.2(a)) or turns out to be
the disagreement above.

---

# TESTS

## What was added, and where

Four `TEST_CASE`s appended to `code/tests/derivation_tests.cpp` — that file only. Nothing under
`code/src`, `code/apps`, `code/CMakeLists.txt`, `derivations/lemmas.md` or `derivations/core.md`
was touched. The file stays standalone doctest.

1. `L1 (lemmas.md) class prediction from the combinatorics vs the numeric class`
2. `L1.2a class-3 harmonics: constant sign on (0, pi], theta = pi included`
3. `L2 (lemmas.md) dim K = 2 rank(D), the (L2.8) generators, and an INDEPENDENT D`
4. `L2 counterexample search: 2000 random (torus, sigma) draws`

Two additions to the file's includes: `method/periodic_jacobian.hpp` and `apps/kill_common.hpp`
(for `reference_cases()` / `make_graph()` / `median_edge_length()`, so the L1 corpus is the shared
kill population rather than a private copy). The K7 torus population, `PState`/`prepare`/`K_at`
and `achievable()` are copied from `code/apps/kill_k7.cpp` / `kill_b3.cpp` so the numbers compare
directly with `results/kill/k7/k7_main.csv`; that copying is marked in the source. The genuinely
independent piece is `dtau_covector`, written from (L2.5) alone.

## Build line

The build line in `derivations/check.md`'s header, plus **`-Icode/apps`** (the only change), and
with the Xcode toolchain compiler by its full path, because a bare `clang++` aborts on this machine
(`xcodebuild` fails to load — the environment defect `lemmas.md` also documents):

```
/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++ \
  -std=c++20 -O2 -arch arm64 \
  -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk \
  -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -Icode/src -Icode/apps \
  code/tests/derivation_tests.cpp code/build/libkiri_core.a \
  -o code/build/derivation_tests && ./code/build/derivation_tests
```

## Whole-file run

```
[doctest] test cases:     37 |     37 passed | 0 failed | 0 skipped
[doctest] assertions: 149115 | 149115 passed | 0 failed |
[doctest] Status: SUCCESS!
```

Baseline before my additions, rebuilt and re-run first so the delta is honest: **33 test cases,
89 074 assertions, 0 failures**. So the four new cases contribute **60 041 assertions**, and no
pre-existing case regressed. Wall time for the whole file: 88 s.

## L1 output

```
L1 corpus: 60 graphs
  L1  triples=1080  candidates=24460884
      classes: zero 223848  class1 23119218  class2 1097596  class3 20222
      coincident-copy pairs: 835272   of which C != 0 (must be 0): 0
        by shared-edge type at v_w:  hinge 447696  split 98352  none 289224
          (of the 'none', f,g share an edge NOT at v_w: 0)
      class predicted from dS == numeric class: 835272 agree, 0 disagree
      NOT coincident but C == 0: 506394  (1e-3 probe: accidental 1120 / structural 505274)
        per-graph persistence (L1.1(iv) criterion): structural 7624  accidental 1116
                                                    rigid/undecidable 20447
      h == 0 pairs: 223848   w IS a or b as an M'-vertex: 223848   not: 0
        index-equal but NOT class zero: 0
      class-3 hinge pairs: 1472   of which <d, x_src - x_dst> = 0: 1472
```

| id | claim | samples | max error | tol | verdict |
|---|---|---|---|---|---|
| L1-a | L1.1(ii) coincident copy ⟹ `h(0) = 0` | 835 272 | 0.000e+00 | 0 | PASS |
| L1-b | L1.2 (L1.8) `p = σ_f⟨d,dS⟩/2` | 835 272 | 2.274e−12 | 1e−9 | PASS |
| L1-c | L1.2 (L1.8) `q = −p` | 835 272 | 3.566e−14 | 1e−9 | PASS |
| L1-d | L1.2 (L1.8) `r = det(d,dS)/2` | 835 272 | 3.411e−13 | 1e−9 | PASS |
| L1-e | L1.2 class from `dS` alone == numeric class | 835 272 | 0.000e+00 | 0 | PASS |
| L1-f | L1.2b type (b) `dS = 2σ_g J(x_src − x_v)` | 223 848 | 1.281e−13 | 1e−9 | PASS |
| L1-g | T1.B/L1.2b type (c) `dS` same for both split endpoints | 98 352 | 7.944e−14 | 1e−9 | PASS |
| L1-h | L1.2b type (b) class 3 ⟺ tested edge ⟂ hinge edge | 1 472 | 0.000e+00 | 0 | PASS |
| L1-i | L1.3 `h ≡ 0` ⟹ `w` is `a` or `b` as an `M′`-vertex | 223 848 | 0.000e+00 | 0 | PASS |
| L1-j | L1.3 `w` is `a` or `b` as an `M′`-vertex ⟹ `h ≡ 0` | 223 848 | 0.000e+00 | 0 | PASS |
| L1-k | H-L1 REFUTED: non-coincident pairs with `C = 0` at every sample | 7 624 | 0.000e+00 | 0 | PASS |
| L1-k2 | L1.1(iii) accidental set (`C = 0` at some samples, not all) | 1 116 | 0.000e+00 | 0 | PASS |
| L1-k3 | L1.1(ii) type (d) as printed: `f,g` share an edge NOT at `v_w` | 835 272 | 0.000e+00 | 0 | PASS |
| L1-l | L1.2a class 3 `= p(1 − cos θ)` identically | 10 000 | 7.767e−16 | 1e−15 | PASS |
| L1-m | L1.2a `sign h = sign p` on `(0,π]`, no zero (`π` included) | 10 000 | 0.000e+00 | 0 | PASS |

**Corpus.** 60 graphs — the 7 non-periodic Phase-2 reference cases from
`kill_common.hpp::reference_cases()` plus `make_graph(id, 18, 46, 220)` for `id = 0 …` until 60 are
usable — `× 18` samples (`X₀` plus Gaussian null-space perturbations at `0.25`), i.e. **1 080
`(graph, σ, t)` triples**, well above the spec's 1 000, and **24 460 884 candidate harmonics** over
every ordered face pair with no broad phase.

**On the spec's "assert zero unexplained instances".** It cannot hold, and `lemmas.md` itself says
why: H-L1 is false (L1.1(iv)). What the test asserts instead:
* the **sound** direction is exact and *is* asserted at zero tolerance — coincident copy ⟹ `C = 0`,
  `835 272 / 835 272`, no exception (`L1-a`);
* the H-L1 remainder is partitioned by the spec's `1e−3` probe **and** by the stronger per-graph
  criterion of L1.1(iv) (`C = 0` at *every* sample of the graph), and both partitions are reported.
  The `1e−3` probe is the weaker of the two and calls far more pairs structural; the per-graph
  criterion is the one whose numbers are comparable with L1.1(iv).

**L1.2a, the `θ = π` question.** 10 000 random exact class-3 harmonics (`q = −p`, `r = 0`), each
sampled at 2 000 angles with `j = 2000` landing **exactly** on `θ = π`: zero sign changes, zero
zeros, `h(π) = 2p` to `1e−14` relative, `harmonic_roots_deflated` empty on `(0, π]`, and
`classify_harmonic` returns `Class3` on all 10 000.

## L2 output

```
    [H2 violated where it matters] squares_3x2            c(Gamma_super)=3   consistency=0.000e+00
    [H2 violated where it matters] squares_3x3            c(Gamma_super)=12  consistency=0.000e+00
    [H2 violated where it matters] snub_square_3x3        c(Gamma_super)=9   consistency=2.750e+01
    [Dind mismatch] snub_square_3x3         k= 18 rankD=2 err=3.715e-02  consistency=2.750e+01
    [H2 violated where it matters] voronoi_torus_1076_n37 c(Gamma_super)=4   consistency=1.078e+01
    [Dind mismatch] voronoi_torus_1076_n37  k= 37 rankD=2 err=1.634e-02  consistency=1.078e+01
  L2  patterns evaluated: 143  (K7 population: 33 of 33)
      rank(A) == 2 rank(D): 143 / 143   failures 0
      dim K in {0,2,4}: failures 0
      rank(D) distribution: 0 -> 8   1 -> 1   2 -> 134
      rank(D) == min(2,k)  [L2.2(e) conjecture]: 142 / 143
      max rel |M_2i   - 2 e_y d_i^T P0^-1| (CIRCULAR, D is read off M_2i) : 2.665e-15
      max rel |M_2i+1 + 2 e_x d_i^T P0^-1| (INDEPENDENT of how D is built): 5.271e-14
      max rel vanishing row of M_j P0                                     : 1.592e-14
      max rel |D from (L2.5) - AchievableSet::D|  (134 patterns)          : 3.715e-02
        restricted to patterns whose period is CONSISTENT                 : 2.664e-13
      consistency > 1e-9 : 2   c(Gamma_super) > 1 : 26   face copies split : 4   D mismatches: 2
      d_tau(1) == 0 (bipartite parity; needed by L2.2(b)/(c))             : 0.000e+00
      path-independence on ker (134 patterns)                             : 1.374e-15
        of those, the two path covectors differ OFF ker: 120
      L2.2(c) rank([L;d_h;d_v]) - rank(L) == rank(D): 108 ok, 0 mismatch
      squares_3x3: k=4 dimK=0 max|D|=3.084e-18 tr K0=0.000000 det K0=1.000000 K0(1,0)=+1.000000
  L2 counterexample search: 2000 drawn, 1999 built;  rank(D) 0/1/2 = 0/0/1999;
                            rank(A) != 2 rank(D): 0;  dim K odd: 0
```

| id | claim | samples | max error | tol | verdict |
|---|---|---|---|---|---|
| L2-a | L2.1 (L2.9) `rank(A) = 2 rank(D)` | 143 | 0.000e+00 | 0 | PASS |
| L2-b | L2.1 `dim 𝒦 ∈ {0,2,4}` — never 1, never 3 | 143 | 0.000e+00 | 0 | PASS |
| L2-c | L2.1 (L2.8) `M_{2i} = 2e_y d_iᵀP₀⁻¹` (**circular**) | 143 | 2.665e−15 | 1e−10 | PASS |
| L2-d | L2.1 (L2.8) `M_{2i+1} = −2e_x d_iᵀP₀⁻¹` (non-circular) | 143 | 5.271e−14 | 1e−10 | PASS |
| L2-e | L2.1 (L2.8) the vanishing row of `M_j P₀` | 143 | 1.592e−14 | 1e−10 | PASS |
| L2-f | L2.0c `D` rebuilt from (L2.5) `==` `AchievableSet::D` (H2 holds) | 132 | 2.664e−13 | 1e−9 | PASS |
| L2-f2 | the `D` mismatches are exactly the inconsistent-period patterns | 143 | 0.000e+00 | 0 | PASS |
| L2-g | L2.2(b)/(c) need `d_τ(1) = 0` (bipartite parity, unprinted) | 134 | 0.000e+00 | 1e−12 | PASS |
| L2-h | L2.0c path-independence holds **on** `ker[L;e_pin]` | 134 | 1.374e−15 | 1e−9 | PASS |
| L2-i | L2.2(c) `rank([L;d_h;d_v]) − rank(L) = rank(D)` | 108 | 0.000e+00 | 0 | PASS |
| L2-j | 2 000-draw counterexample search for `rank(A) = 2rank(D)` | 1 999 | 0.000e+00 | 0 | PASS |

**Corpus.** The **33 K7 patterns** (`l2_k7_population()`, byte-identical to `kill_k7.cpp`'s
`population()`; all 33 build and solve) **plus 110 fresh random Voronoi tori** of 12–45 sites with
`quotient_sigma` max-cut σ, of which 110 build → **143 patterns**. Ranks are dense `JacobiSVD` with
threshold `1e−10·σ_max` (floor `1e−13`), as the spec asks. The covector test (L2.10) runs on the
**108** patterns with `n_q ≤ 80`; on the larger tori it is not measured, and that is a limitation,
not a pass.

**Counterexample search.** 2 000 random `(torus, σ)` draws (6–16 sites, fresh σ per draw), 1 999
built; `rank(A) = 2 rank(D)` on **1 999 / 1 999**, `dim 𝒦` even on all. **No counterexample.**
All 1 999 have `rank(D) = 2`, which is further evidence for — not a proof of — L2.2(e).

**`rank(D)` distribution over the 143.** `0 → 8`, `1 → 1`, `2 → 134`. The eight zeros are the seven
`k = 0` patterns (`squares_2x2`, all three `triangles`, all three `kagome`) plus `squares_3x3` with
`k = 4`; the single `1` is `squares_3x2` with `k = 1`. This reproduces `lemmas.md`'s Check L2
distribution (`0 → 8`, `1 → 1`, `2 → 124`) with 10 more tori, and reproduces L2.2(e)'s 132/133 as
**142/143 with the same sole exception**.

**What I could not test.**
* The near-tolerance population of the Deriver's weakest step 1 (a structural `C ≡ 0` pair landing
  within a decade of `1e−11·scale²`, where `validity_certificate` could return `noroot = true`
  unsoundly). My tests use the same `1e−11·scale²` rule and therefore inherit the question rather
  than answering it. This is the **only** open soundness question I am leaving.
* Whether a pattern realising `K = αI` with `α < 0` exists (D-L2-6). The proof defect stands either
  way; the question is whether the case is reachable.
* `rank(D)` via (L2.10) on the 35 patterns with `n_q > 80`.
* Every `[N]` number inside `lemmas.md`'s own Check L1 / Check L2 sections: I ran my own programs on
  my own corpora, so my counts are comparable, not identical, except where noted (the accidental set
  reproduces exactly at **1 116**).

---

# OVERALL VERDICT

**Lemma L1.1 — VERIFIED.** (i) and (ii)'s implication are proved and I proved them blind; (iii) is
proved and in fact needs no hypothesis; (iv) is an honest measurement of H-L1's failure and I
reproduce its accidental set exactly (1 116). Two defects, neither mathematical: the type-(d) row
should read "no edge **containing `v_w`**" (D-L1-1, 0 occurrences measured), and (iv)'s counts do
not state whether they are pairs or `(pair, sample)` instances — on the reading that makes two of
the three numbers agree with my run, the third is off by a factor of about 20 (D-L1-4). The lemma's
load-bearing content, that `h(0) ≡ 0` is structural for coincident-copy pairs and that a purely
combinatorial classifier under-counts the `C = 0` class in the *conservative* direction, is correct
and is the right thing for the certificate to rely on.

**Lemma L1.2 — VERIFIED.** The one lemma I would sign without reservation. (L1.8) is right in both
endpoint cases — I derived it independently before reading the proof, and the `v = v_b` case, where
the cross terms must cancel, is the place a sign could have hidden and does not. It cross-checks
exactly against `core.md` (T5.2)/(T5.4) for split edges and against `zero_plus.hpp`'s `dS_e`. The
trichotomy is exhaustive; (L1.2a) reproduces Lemma T5.1e and my attack on its "never a contact"
reading failed for the right reason (a permanent touching is `h ≡ 0`, which is a *different* class);
and step 5 is a real contribution, turning `check.md` R3.1a's hand-added `A ≠ 0` side-condition into
the consequence `dS ≠ 0`. Measured `835 272 / 835 272` with zero class disagreements. Only (L1.2c)'s
"codimension ≥ 1" is unquantified, and (L1.2c) is not load-bearing.

**Lemma L1.3 — VERIFIED-UNDER-HYPOTHESIS**, the hypothesis being **"every split cut opens"
(`Δu ≠ 0` on every split-adjacent face pair), plus the absence of the type-(d) coincidence**
`2Δu = (σ_g − σ_f)x_{v_w}`. The equivalence `1 ⇔ 2 ⇔ 3` is fully proved, blind and again in the
Deriver's frame, and the hinge specialization (L1.9) is right. What is **not** proved is the boxed
word "**exactly**": the Deriver's own proof steps 5 and 6 concede that split-adjacent pairs with
`Δu = 0`, and non-adjacent pairs at isolated `t`, also satisfy condition 3 — and step 5's reply
(such designs are rejected by the certificate) is about the certificate, not about the pair. The
named hypothesis is not free: `Δu = 0` is exactly the non-opening split cut that
`zero_plus.hpp`/K5/F30 report as the *observed* failure mode. Relatedly, the box's "decided by
`(M, σ)` … not by `t`" is false for the split branch, `Δu` being linear in `X`. On the corpus the
characterisation is exact in **both** directions (`223 848 / 223 848`, tests `L1-i`/`L1-j`), and it
is sharper than the lemma prints: the class is decided by an `M′`-vertex **index comparison**, which
is what `contact.cpp` already does. So the engineering claim is safe and the theorem needs its
hypothesis printed. **DISPUTED on the word "exactly" (D-L1-3, D-L1-2).**

**Corollary L1.4 — VERIFIED**, inheriting L1.1 and L1.3's caveats. Its (iii) is the most honest
paragraph in the file and its conservativeness argument is correct as an argument about the *exact*
predicate. The Deriver's own flagged residue — the near-threshold numeric population — is real, is
not closed by his tests or by mine, and is the one place I would look next for an unsoundness.

**Lemma L2.1 — VERIFIED.** Both (L2.8) generators, both signs, and (L2.9) are proved; I re-derived
all of it blind and reached the same formulas. The spec's attack on the factorisation fails for a
reason worth recording: the `x`-copy and the `y`-copy of one null vector occupy **complementary
rows** of `M_j P₀`, so no `P₀` — square lattice or otherwise — can collapse them, and the rank is
forced even. The measurement now has a non-circular leg it did not have: `D` rebuilt directly from
(L2.5) matches `AchievableSet::D` to `2.66e−13` on 132 of 134 patterns, which closes the Deriver's
own "weakest step 3". **VERIFIED-UNDER-HYPOTHESIS in one respect (D-L2-5):** H2 must be stated for
the graph the period is *measured* on, not only for the quotient. On `snub_square_3x3` — a pattern
inside Check L2's own "133/133" — the `3×3` super patch's hinge graph has 9 components,
`PeriodicJac::consistency = 27.5`, and the code's `D` disagrees with (L2.5)'s by `3.7e−2`. Check L2
should read **132/133 measured + 1 excluded**, and `achievable()` should be gated on `consistency`.

**Lemma L2.2 — DISPUTED (D-L2-2), conclusions VERIFIED.** (a) is trivial by shape and says so.
(b) and (c) are **true** and their proofs contain a **false step**: `{φ_i}` is treated as a basis of
`ker L` when it is a basis of `ker[L; e_pin]`, one dimension smaller, because the constant vector
lies in `ker L` and is removed by the pin. The conclusions survive only via `d_τ(1) = 0`, which the
file never states; it follows in two lines from the bipartiteness of `Γ` (every lifted walk
`f₀ → f₀+τ` has even length) and I measure it exactly `0.000e+00` on 134 patterns. (c)'s duality
argument is otherwise correct and I reproduce (L2.10) on 108/108 patterns. (d) is a correct and
valuable worked case whose stated evidence is one number short of its conclusion (`tr` and `det` do
not separate `J` from `−J`; the sign is `+1`, so the claim is true). (e) is **ILL-POSED as a
statement** — its exception clause is not a predicate — but is correctly labelled a conjecture,
correctly unused, and its evidence reproduces at 142/143 with the same exception.

**Corollary L2.3 — DISPUTED (D-L2-6), statement VERIFIED.** Both implications and (L2.11) are
right. The non-singularity paragraph is not: it concedes in its last clause exactly what it sets out
to exclude, and `K = −I` gives `J(π/2) = 0` inside H6's range `(0, π]`. Either read the conclusion
in the linear span `span{I,J}`, which is all L2.3 is used for, or add the side-condition
`K ∉ {αI : α < 0}`.

**Nothing is REFUTED.** No counterexample to any statement was found in either pass, in 24 460 884
candidate harmonics, 143 periodic patterns, or 2 000 fresh `(torus, σ)` draws.
