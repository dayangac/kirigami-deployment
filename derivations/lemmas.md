# lemmas.md — two results that `core.md` / `KILL_REPORT.md` currently state as MEASURED

Deriver-L, mission 2 / WP1. Companion to `derivations/core.md`; nothing here edits that file,
`code/src` or `code/tests`. Notation is `core.md` §0 throughout (`c = cos(θ/2)`, `s = sin(θ/2)`,
`τ = tan(θ/2)`, `J` the `+π/2` rotation, `σ : F → {±1}`, `u` the face potential of T1 Step 3,
`Y_θ = cC + sS` with `C_{(v,f)} = x_v`, `S_{(v,f)} = J(2u_f − σ_f x_v)`).

Tags, as in the technical report: **[D]** definition · **[F]** cited fact, with its source ·
**[A]** algebra a reader can verify by hand · **[N]** numerical measurement. Nothing is
"obvious"; every step is one of those four.

## Standing hypotheses

Imported from `core.md` by name, unchanged:

* **H1** `M = (V,E,F)` is a straight-line embedded planar graph, `σ` given, faces stored
  counter-clockwise (core.md 0.1–0.2).
* **H2** `Γ` is connected, `c(Γ) = 1` (core.md T1 hypotheses; `[F, F22]` 50/50 under Eq. (1)).
* **H3** `X ∈ 𝕏 = X₀ + span Φ`, i.e. Eq. (2) holds on every hole preimage, so the face potential
  `u` exists and is `θ`-free (core.md T1 Steps 3–6).
* **H4** every face polygon of `M` is simple and positively oriented at `θ = 0` (core.md T4
  standing hypotheses; by T1.D this then holds at every `θ`).
* **H5** no edge of `M` has zero length: `x_{dst} ≠ x_{src}` for every edge.
* **H6** the deployment runs over `θ ∈ (0, π]` (core.md deviation 9).

New, and used only where named:

* **H-L1 (flat genericity).** For every candidate pair `π = (w,(a,b))` whose three **source**
  vertices `ρ(a), ρ(b), ρ(w)` are three *distinct* vertices of `M`, the quadratic polynomial
  `t ↦ det( x_{ρ(b)} − x_{ρ(a)} , x_{ρ(w)} − x_{ρ(a)} )` on the shape space is **not** identically
  zero. **H-L1 is FALSE on the measured corpus** — see L1.1(iv) and Check L1 — and only the
  *converse half* of Lemma L1.1 needs it. Everything the certificate actually uses (L1.1(i)–(ii),
  L1.2, L1.3, L1.4) is proved without it.
* **H-L3 (every split cut opens, and no accidental welding) — [R2 correction: D-L1-3], new in
  round 2.** For every pair of faces `f, g` adjacent across a **split** edge, `Δu = u_g − u_f ≠ 0`;
  and for every pair of faces meeting only at a vertex `v_w` (type (d) of L1.1(ii)),
  `2Δu ≠ (σ_g − σ_f) x_{v_w}`. Used **only** to close Lemma L1.3's catalogue (L1.3-cor) and the
  first bullet of L1.4(i). Without it the catalogue is still correct in the direction the
  certificate uses, but not exhaustive: split-adjacent pairs with a non-opening cut, and type-(d)
  pairs, also have `h ≡ 0`. The first clause is exactly the `0⁺` separation predicate of
  `zero_plus.hpp`; `[F, K5 / F30]` reports non-opening split cuts as the observed failure mode, so
  H-L3 fails precisely on designs with `Θ_max = 0`, where no contact calculus is needed.
* **H-L2 (non-degenerate lattice).** `P₀ = [t_h t_v]` is invertible. This is automatic for a
  lattice: two independent periods. Without it `K = Q P₀⁻¹` is undefined and `J(θ)` is not a
  Jacobian at all.

---

# Part L1 — the degeneracy catalogue of the contact harmonics

## L1.0 The shared computation: one orthogonal frame kills all four cases

**[D] L1.0a (candidate pair).** A *candidate pair* is `π = (w, (a,b))` where `(a,b)` is an
`M′`-edge of some face `f` (two cyclically consecutive corners of `f`) and `w` is a corner of a
face `g ≠ f`. Write `ρ : V′ → V` for the source map (`cut.hpp::prime_to_original`) and
`v_a = ρ(a)`, `v_b = ρ(b)`, `v_w = ρ(w)`. The *orientation harmonic* of `π` is
`h_π(θ) = det(Y_b − Y_a, Y_w − Y_a) = p + q cos θ + r sin θ`, coefficients by (T3.2).

**[A] L1.0b (the frame).** Because `a` and `b` are corners of the **same** face `f`,
core.md T3.3's computation gives, with `d := x_{v_b} − x_{v_a} ≠ 0` (H5),

```
        C_ab = d ,          S_ab = − σ_f J d ,                                       (L1.1)
```

so `{C_ab, S_ab}` is an **orthogonal basis of `R²`** with

```
        D₀ := det(C_ab, S_ab) = − σ_f det(d, J d) = − σ_f |d|²  ≠ 0 .                (L1.2)
```

Expand the other two columns in that basis:

```
        C_aw = α C_ab + β S_ab ,        S_aw = γ C_ab + δ S_ab .                     (L1.3)
```

**[A] L1.0c (the four coefficients in the frame).** Substituting (L1.3) into (T3.2) and using
`det(C_ab,C_ab) = det(S_ab,S_ab) = 0`, `det(S_ab,C_ab) = −D₀`:

```
   ┌───────────────────────────────────────────────────────────────────────────────┐
   │  C := p + q = det(C_ab, C_aw) =  β D₀                                         │
   │  A := p − q = det(S_ab, S_aw) = −γ D₀                                         │
   │  B := 2r    = det(C_ab,S_aw) + det(S_ab,C_aw) = (δ − α) D₀        (L1.4)      │
   └───────────────────────────────────────────────────────────────────────────────┘
```

Since `D₀ ≠ 0`, the three-class split of core.md T5.2b.2 and the identity test of T3.H.1 become
**four linear conditions on `(α, β, γ, δ)`**:

```
   class 1   β ≠ 0
   class 2   β = 0 ,  δ ≠ α
   class 3   β = 0 ,  δ = α ,  γ ≠ 0                                                (L1.5)
   zero      β = 0 ,  δ = α ,  γ = 0
```

**[A] L1.0d (chart-free class definition).** `C = h(0)`, `B = 2h′(0)` and, when `C = B = 0`,
`A = 2h″(0)`: indeed `h(0) = p+q`, `h′(0) = r`, `h″(0) = −q = (A−C)/2`. So the class split is a
statement about the **order of vanishing of `h` at `θ = 0`** and never mentions `τ = tan(θ/2)`.
This removes the `θ = π` singularity of the `τ` chart (`check.md` D10) from the *classification*;
it still lives in the *root formulas*, where core.md already handles it.

## Lemma L1.1 — structural `h(0) = p + q = 0`

**Statement.**

**(i)** For every candidate pair and every `X`,

```
        h_π(0)  =  p + q  =  det( x_{v_b} − x_{v_a} ,  x_{v_w} − x_{v_a} ) ,        (L1.6)
```

the orientation determinant of the three **source** vertices in the flat state. In particular
`h_π(0)` depends on the pair only through `(v_a, v_b, v_w)`, not through which `M′`-copies were
taken, and it is a quadratic polynomial in the design coordinate `t`.

**(ii)** *(sufficient, purely combinatorial)* If `v_w = v_a` or `v_w = v_b` — the pair is a
**coincident-copy pair** — then `h_π(0) = 0` for **every** `X ∈ R^{N×2}`, a fortiori on all of
`𝕏`. The coincident-copy pairs fall into exactly four combinatorial types:

```
   (a) f, g hinge-adjacent across e ,  v_w = src(e)      permanent incidence (Lemma L1.3)
   (b) f, g hinge-adjacent across e ,  v_w = dst(e)      the duplicated hinge endpoint
   (c) f, g split-adjacent across e ,  v_w ∈ {a_e, b_e}  the split-duplicate pair of T5.3
   (d) f, g share the vertex v_w but no edge CONTAINING v_w   copies around a cut cone
```

*[R2 correction: D-L1-1]* Row (d) previously read "share the vertex `v_w` but no **edge**", which
left uncovered the configuration "`f, g` adjacent across an edge `e`, and additionally meeting at
some `v_w ∉ e`". Proof step 3 always had the right dichotomy; only the table row was wrong. With
"no edge containing `v_w`" the four types are exhaustive, and they are **not** disjoint (a pair can
be of type (c) at both endpoints at once). Measured 0 occurrences of the uncovered configuration
(Checker test `L1-k3`, 835 272 coincident-copy pairs).

**(iii)** *(converse — unconditional; H-L1 is NOT used here)* **[R2 correction: D-L1-6]**
If `v_a, v_b, v_w` are three distinct vertices, then either
`h_π(0)` vanishes identically on `𝕏` — in which case the three flat points are collinear for
**every** design — or `{ t : h_π(0; t) = 0 }` is the zero set of a non-zero polynomial of degree
`≤ 2`, hence a proper algebraic subset of `R^{2k}` of Lebesgue measure zero: an **accidental**
degeneracy. The dichotomy itself needs no hypothesis. **H-L1 is the separate assertion that the
first alternative never occurs**, and it is that assertion — not (iii) — that (iv) refutes. Round 1
printed (iii) as "(converse, under H-L1)", which mislabelled an unconditional statement.

**(iv)** *(H-L1 fails, and how)* **[R2 correction: D-L1-4]** Every count in this paragraph and in
Check L1 is an **instance** count — one unit per `(candidate pair, shape-space sample)` — not a
count of distinct pairs; `check_l1.cpp` increments once per sample. The Checker's independent run
reports **pair** counts, which is why its structural figure (7 624 pairs) is ~20× smaller than the
174 354 instances below while its accidental figure (1 116) coincides with mine, those pairs each
occurring in a single sample. On the corpus of Check L1 (84 graphs, 30 119 368 candidate-harmonic
**instances**) **197 721** instances (0.66 %) have `|h_π(0)| ≤ tol` with `v_w ∉ {v_a, v_b}`. Of those,
**174 354** vanish at **every** shape-space sample of their graph, so they are structural, not
accidental: **135 410** have all three source vertices *frozen* (`Φ`-row zero — collinear triples
pinned by the fixed boundary), and **38 944** have one or two free vertices that the shape space
nevertheless keeps collinear (every witness the program printed is in `tiling_3_4_3_12`; I did not
enumerate the rest, so "all in `tiling_3_4_3_12`" is not claimed). A further **22 251** sit on graphs with
`dim_null = 0`, where there is no shape space to test. Only **1 116** of 30 119 368 (`3.7e−5`)
vanish at some samples and not others — the accidental set of (iii). The rate quoted for
"**H-L1 is false**" is the *structural* one, `174 354 / 30 119 368 = 0.579 %`; the 0.66 % headline
above additionally includes the accidental and the `dim_null = 0` instances. **[R2 correction:
D-L1-4]**

**Proof.**

1. **[A]** `Y_θ(0) = C` by (T1.5), and `C_{(v,f)} = x_v` for every copy. Hence
   `h_π(0) = det(C_b − C_a, C_w − C_a) = det(x_{v_b} − x_{v_a}, x_{v_w} − x_{v_a})`, which is
   (L1.6). Independently, evaluating `p + q cos θ + r sin θ` at `θ = 0` gives `p + q`. ∎(i)
2. **[A]** If `v_w = v_a` the second argument of the determinant is `0`; if `v_w = v_b` the two
   arguments are equal. Either way the determinant vanishes for arbitrary vertex positions, so no
   hypothesis on `X` is used. ∎(ii, first half)
3. **[F, core.md 0.2 / F11]** The four types are exhaustive: `w` is a corner of `g` at `v_w` and
   one of `a, b` is a corner of `f` at the same `v_w`, so `f` and `g` are two faces meeting at the
   vertex `v_w`. Two faces of a planar graph meeting at a vertex either share an edge at that
   vertex — hinge (`σ_f = −σ_g`) or split (`σ_f = σ_g`) by 0.2, and in the hinge case `v_w` is
   `src(e)` or `dst(e)` since it is an endpoint of `e` — or share **no edge containing `v_w`**,
   which is type (d) as amended. ∎(ii)
4. **[A]** By (i), `h_π(0)` is `det` of two differences of rows of `X = X₀ + Φ t`, hence a
   polynomial of degree `≤ 2` in `t`. A non-zero real polynomial vanishes on a set of Lebesgue
   measure zero (a proper algebraic subset). ∎(iii)
5. **[N]** Check L1, the block "L1.1". ∎(iv)

**What it replaces.** core.md T3.H.5's sentence "`g(0) = p + q = 0` … says the triple `(a,b,w)` is
collinear in the *flat* state, which holds *identically on all of `𝕏`* whenever `w` and `a` (or `w`
and `b`) are two `M′`-copies of the same vertex of `M`", and T5.2b.2's "Membership in class 2 or
class 3 is decided by the *combinatorics* `(M,σ)`". Both are correct in the direction they are
used, and (iv) says the *converse* is not: the combinatorial rule is a strict subset of the
structural `C ≡ 0` set, so a classifier that reads the combinatorics alone must fall back on the
numeric test for the remaining 0.6 %. core.md's own caveat (1) in T5.2b.2 anticipates this and is
now quantified.

## Lemma L1.2 — the class of a coincident-copy pair, and why class 3 is never a contact

**Statement.** Let `π = (w,(a,b))` be a coincident-copy pair (L1.1(ii)) with common source vertex
`v` (`v = v_a` or `v = v_b`), let `a*` be the copy of `v` inside `f` (so `a* = a` or `a* = b`), and
put

```
        d  :=  x_{v_b} − x_{v_a}  ≠ 0 ,        dS  :=  S_w − S_{a*}  ∈ R² .          (L1.7)
```

`dS` is `2 ∂/∂θ|₀` of the separation of the two copies of `v`: `Y_w(θ) − Y_{a*}(θ) = sin(θ/2)·dS`
exactly, for all `θ`. Then

```
   ┌───────────────────────────────────────────────────────────────────────────────┐
   │   p = ½ σ_f ⟨d, dS⟩ ,     q = −p ,     r = ½ det(d, dS) ,                      │
   │   so   C = 0 ,   A = σ_f ⟨d, dS⟩ ,   B = det(d, dS) .              (L1.8)      │
   └───────────────────────────────────────────────────────────────────────────────┘
```

Consequently the whole class of `π` is decided by the single vector `dS`:

```
   dS = 0                                →  h_π ≡ 0          (Lemma L1.3)
   dS ≠ 0 and dS ∥ d                     →  class 3
   det(d, dS) ≠ 0                        →  class 2, with the one further root
                                            τ* = −B/A = − det(d,dS) / ( σ_f ⟨d,dS⟩ )
```

and:

**(L1.2a) class 3 is never a contact.** If `h(0) = h′(0) = 0` and `h″(0) ≠ 0` then
`h(θ) = p(1 − cos θ)` with `p = h″(0) ≠ 0`, which has **constant sign `sign p` on `(0, 2π)`** and
its only zero in `[0, π]` is `θ = 0`. Hence `π` contributes no candidate contact angle whatsoever,
and the deflated atom for it is the constant `true`.

**(L1.2b) what class 3 is, geometrically.** `dS ∥ d` says the two copies of `v` separate, to first
order and in fact at every `θ`, **along the direction of the tested edge**: the moving copy slides
along the line of `(a,b)` instead of crossing it. Explicitly:

```
   type (b), hinge, v = dst(e):   dS = 2 σ_g J ( x_{src(e)} − x_v )
                                  class 3  ⟺  ⟨ d , x_{src(e)} − x_{dst(e)} ⟩ = 0
                                            i.e. the tested edge ⟂ the hinge edge
   type (c), split:               dS = 2 J Δu       (T1.8; the same vector for both endpoints)
                                  class 3  ⟺  ⟨ d , Δu ⟩ = 0 ,  which is exactly T5.3's r = 0
                                            (core.md T5.H.4: "the flat root is double")
   type (d), vertex only:         dS = J ( 2Δu − (σ_g − σ_f) x_v ) ,   Δu = u_g − u_f
```

**(L1.2c) how many.** Class 3 is **one** polynomial equation, `det(d, dS) = 0`, imposed on top of
the coincident-copy condition; on a coincident-copy pair `det(d,dS)` is a quadratic polynomial in
`t`, so class 3 is a codimension-`≥ 1` condition — *accidental* for a generic design, but *forced*
whenever the pattern's geometry makes the two directions orthogonal, which is what the reference
tilings do. This is why the count is non-zero, `ε`-independent, and small (core.md's 1 138 of
2 145 387).

**Proof.**

1. **[A]** `v_w = v` and `ρ(a*) = v`, so `C_w = C_{a*} = x_v` and `Y_w − Y_{a*} = s (S_w − S_{a*})`
   by (T1.5): the `cos(θ/2)` part cancels identically. That is the sentence after (L1.7), and it
   is the same computation as `zero_plus.hpp`'s header (`dC_e = 0`, `dS_e`) generalised from split
   edges to every coincident copy.
2. **[A]** In the frame of L1.0b, `C_aw = C_w − C_a`. If `v = v_a` then `C_aw = 0`, i.e.
   `α = β = 0`, and `S_aw = dS`. If `v = v_b` then `C_aw = C_ab`, i.e. `α = 1, β = 0`, and
   `S_aw = S_ab + dS`. In both cases `β = 0`, so `C = 0` by (L1.4) — re-deriving L1.1(ii) — and, in
   both cases, writing `dS = γ C_ab + δ' S_ab`, one gets `A = −γ D₀` and `B = δ' D₀`.
3. **[A]** `det(d, dS) = det(C_ab, γ C_ab + δ' S_ab) = δ' D₀ = B = 2r`, giving `r = ½ det(d,dS)`.
4. **[A]** `⟨d, dS⟩ = ⟨C_ab, γ C_ab + δ' S_ab⟩ = γ |d|²` because `S_ab = −σ_f J d ⟂ d`; and
   `|d|² = −σ_f D₀` by (L1.2), so `⟨d,dS⟩ = −σ_f γ D₀ = σ_f A = 2 σ_f p`, i.e.
   `p = ½ σ_f ⟨d,dS⟩`. With `C = p + q = 0` this gives `q = −p`. That is (L1.8).
5. **[A]** From (L1.8): `A = B = 0 ⟺ ⟨d,dS⟩ = det(d,dS) = 0 ⟺ dS = 0` (two independent linear
   functionals of `dS`, since `d ≠ 0`); `B = 0` with `dS ≠ 0` `⟺ dS ∥ d`, and then
   `⟨d,dS⟩ ≠ 0` automatically, so `A ≠ 0` and the side condition of (T5.1e′) holds
   **for free** — the case `A = C = B = 0` that `check.md` R3.1a had to exclude by hand is here
   exactly `dS = 0`. ∎ (the trichotomy)
6. **[A]** For class 2, `g(τ) = Bτ + Aτ²` (core.md T5.2b.2), whose non-zero root is
   `τ* = −B/A = − det(d,dS) / (σ_f ⟨d,dS⟩)`. For a split edge this is core.md (T5.4)
   `τ* = −r/p` verbatim, since there `p = −σ_f det(d,Δu)` and `r = ⟨d,Δu⟩`.
7. **[A] (L1.2a)** `h(0) = 0` and `h′(0) = 0` give `q = −p`, `r = 0`, so
   `h(θ) = p − p cos θ = p(1 − cos θ)` and `h″(0) = −q = p`. For `θ ∈ (0, 2π)`, `cos θ < 1`
   because `cos θ = 1` exactly on `2πZ`; hence `1 − cos θ > 0` and `h(θ)` has the sign of `p ≠ 0`
   throughout, with no zero. The argument uses only the cosine, never `τ`, so nothing about it
   degenerates at `θ = π` (`check.md` D10). ∎(L1.2a)
8. **[A] (L1.2b)** For type (b), `f` and `g` are hinge-adjacent across `e`, so
   `Δu = u_g − u_f = σ_g x_{src(e)}` by (T1.3) and `σ_g − σ_f = 2σ_g`. Then
   `dS = S_w − S_{a*} = J(2Δu − (σ_g − σ_f) x_v) = 2σ_g J(x_{src(e)} − x_v)`, and with
   `v = dst(e)` this is `2σ_g J(x_{src} − x_{dst})`. `det(d, J z) = ⟨d, z⟩`, so
   `B = det(d,dS) = 2σ_g ⟨d, x_{src} − x_{dst}⟩`. For type (c), `σ_f = σ_g` so the second term
   drops and `dS = 2 J Δu`, which is (T1.8) and is independent of which endpoint of the split edge
   `v` is — core.md T1.B. Then `B = 2 det(d, JΔu) = 2⟨d, Δu⟩ = 2r`, matching (T5.2). ∎(L1.2b)
9. **[N]** (L1.8), the trichotomy and both closed forms are measured in Check L1: `801 816`
   coincident-copy pairs, **0** class disagreements, `max |p − ½σ_f⟨d,dS⟩| = 1.99e−12`,
   `max |r − ½det(d,dS)| = 3.41e−13`, hinge `dS` formula `1.00e−13` over `479 784` pairs, split
   `dS` endpoint-independence `8.64e−14` over `159 816` pairs, and `160 / 160` of the hinge-type
   class-3 pairs satisfy the perpendicularity `⟨d, x_src − x_v⟩ = 0`. ∎

**What it replaces.** core.md Lemma T5.1e (which is correct, and is step 7 here) plus the two
sentences that surround it: "class 3 is the sub-case in which the flat configuration additionally
satisfies `r ≡ 0`, which by (T3.2) is again a statement about which `M′`-copies coincide, not about
`t`" (T5.2b.2) — that sentence is **wrong as printed**: `r = ½ det(d, dS)` *is* a function of `t`,
and Check L1 finds 152 class-3 occurrences that appear at some `t` and not others. What is true, and
is proved here, is that the *vanishing of `C`* is combinatorial while the *class-3 refinement* is
one extra polynomial equation. Also replaces the side-condition discussion of (T5.1e′) /
`check.md` R3.1a: `A ≠ 0` is not an extra hypothesis, it is `dS ≠ 0` (step 5).

## Lemma L1.3 — the identically-zero harmonics are a combinatorial class

**Statement.** For a candidate pair `π = (w,(a,b))` the following are equivalent:

1. `h_π ≡ 0` on `R` (equivalently `p = q = r = 0`);
2. there is a **constant** `λ ∈ R` with `Y_w(θ) − Y_a(θ) = λ ( Y_b(θ) − Y_a(θ) )` for every `θ`;
3. `x_{v_w} = x_{v_a} + λ ( x_{v_b} − x_{v_a} )` **and** `2 Δu = (σ_g − σ_f) x_{v_w}`, where
   `Δu = u_g − u_f`.

If moreover `f` and `g` are **hinge-adjacent** across `e`, condition 3 reads

```
        x_{v_w} = x_{src(e)}      and      x_{src(e)} ∈ the line of (x_{v_a}, x_{v_b}) ;   (L1.9)
```

if `f` and `g` are **split-adjacent**, condition 3 reads **[R2 correction: D-L1-2]**

```
        Δu = 0        and       x_{v_w} ∈ the line of (x_{v_a}, x_{v_b}) ;             (L1.10)
```

`Δu = 0` says the two faces carry the same face potential, so by (T1.8) the split cut never opens.
Round 1 printed only the first half of (L1.10); the collinearity half was in proof step 5 but not in
the statement.

**[R2 correction: D-L1-3] — the round-1 box, and what replaces it.** Round 1 closed the statement
with "*so the identically-zero pairs are **exactly** the hinge-point pairs*", and called the
description "combinatorial … not by `t`". **Both claims are withdrawn.** The equivalence
`1 ⇔ 2 ⇔ 3` is unconditional and stands; the word "exactly" does not, because proof steps 5 and 6
themselves exhibit two further families satisfying condition 3 — split-adjacent pairs with `Δu = 0`
plus collinearity (L1.10), and non-adjacent pairs with `2Δu = (σ_g − σ_f) x_{v_w}` plus
collinearity. Round 1's reply to the first family ("such a design is rejected by the very first
predicate of the certificate") is a statement about the *certificate*, not about whether the pair
has `h ≡ 0`, so it does not support "exactly". And `Δu = u_g − u_f` is **linear in `X`** (T1 Step 8),
so `Δu = 0` is a condition on `t`, not a combinatorial one; "not by `t`" was false for the split
branch. What can be stated is:

**(L1.3-cor) [A] Under H-L3 the catalogue closes.** Add the hypothesis

* **H-L3 (every split cut opens, and no accidental welding).** `Δu = u_g − u_f ≠ 0` for every pair
  of faces `f, g` adjacent across a **split** edge, and `2Δu ≠ (σ_g − σ_f) x_{v_w}` for every pair
  of faces meeting only at a vertex `v_w` (type (d)).

Then the identically-zero pairs are exactly

> the pairs in which `w` is the copy in `g` of the **hinge point** `src(e)` of a hinge edge `e`
> shared by `f` and `g`, and that hinge point lies on the line of the tested edge `(a,b)` of `f` —
> in particular every pair of type (a) of L1.1(ii), where the hinge point *is* an endpoint of
> `(a,b)`, so the collinearity holds with `λ ∈ {0,1}` for free.

*Proof.* By 1 ⇔ 3, an identically-zero pair satisfies condition 3. Its two faces meet at `v_w`, so
by L1.1(ii) as amended they share a hinge edge at `v_w`, share a split edge at `v_w`, or share no
edge containing `v_w`. H-L3's first clause kills the split case via (L1.10), its second clause kills
type (d), and the hinge case is (L1.9). ∎

**H-L3 is not free.** Its first clause is exactly the predicate `zero_plus.hpp` is written to
measure (`q_e ∝ ⟨d, Δu⟩`-type first-order separation), and `[F, K5 / F30 / KILL_REPORT]` reports
non-opening split cuts as the *observed* failure mode of projected designs. So H-L3 holds on every
design that survives the `0⁺` test and fails on exactly the designs that do not — which are the ones
with `Θ_max = 0`, where no contact calculus is needed anyway. Without H-L3 the catalogue is still
correct in the direction the certificate uses (`h ≡ 0 ⟹` condition 3), just not exhaustive.

**What is measured, without any hypothesis.** On the Check L1 corpus the catalogue is exact in both
directions: `239 892 / 239 892` identically-zero instances are hinge-`src` pairs, `0` are
split-adjacent, `0` are non-adjacent — and, sharper, all of them have `w == a` or `w == b` **as
`M′`-vertices** (an index comparison). The Checker reproduces both directions on its own corpus,
`223 848 / 223 848`.

**Proof.**

1. **[A]** `(1) ⟹ (2)`. In the frame of L1.0b, `p = q = r = 0` is `C = A = B = 0`, i.e. by (L1.4)
   `β = 0`, `γ = 0`, `δ = α`. So `C_aw = α C_ab` and `S_aw = α S_ab` with the **same** `α`.
   Since `Y_w − Y_a = c C_aw + s S_aw` and `Y_b − Y_a = c C_ab + s S_ab`, this is exactly (2) with
   `λ = α`. (This is where the *orthogonality* of the frame does the work: without it, `(1)`
   would only give pointwise proportionality with a possibly `θ`-dependent factor.)
2. **[A]** `(2) ⟹ (1)`. `det(λE, E) = 0` for every `θ`.
3. **[A]** `(2) ⟺ (3)`. `C_aw = λ C_ab` is `x_{v_w} − x_{v_a} = λ(x_{v_b} − x_{v_a})`. For the
   `S` part, `S_aw = J(2u_g − σ_g x_{v_w}) − J(2u_f − σ_f x_{v_a}) = J(2Δu − σ_g x_{v_w} + σ_f x_{v_a})`
   and `λ S_ab = −λ σ_f J(x_{v_b} − x_{v_a}) = −σ_f J(x_{v_w} − x_{v_a})`. Equating and cancelling
   the invertible `J`: `2Δu − σ_g x_{v_w} + σ_f x_{v_a} = −σ_f x_{v_w} + σ_f x_{v_a}`, i.e.
   `2Δu = (σ_g − σ_f) x_{v_w}`.
4. **[A]** *Hinge-adjacent.* `(T1.3)` gives `Δu = σ_g x_{src(e)}` and `σ_g − σ_f = 2σ_g`, so
   condition 3's second half is `2σ_g x_{src(e)} = 2σ_g x_{v_w}`, i.e. `x_{v_w} = x_{src(e)}`; the
   first half is the collinearity. When `v_w = src(e)` is itself an endpoint of the tested edge
   `(a,b)` — type (a) — the collinearity holds with `λ ∈ {0,1}` for free.
5. **[A]** *Split-adjacent.* `σ_g = σ_f`, so condition 3 becomes (L1.10): `Δu = 0` **and** the
   collinearity. By (T1.8) the offset between the two duplicates is `2 sin(θ/2) J Δu ≡ 0`, so the
   two copies of the whole split edge coincide at every `θ` and `Θ_max = 0` by inspection.
   **[R2 correction: D-L1-3]** Round 1 added "such a design is rejected by the very first predicate
   of the certificate" and read that as an exclusion. It is not: it is a statement about the
   certificate, not about whether the pair has `h ≡ 0`. This family is excluded only by H-L3's
   first clause, and is measured empty on the corpus.
6. **[A]** *Non-adjacent (type (d)).* `Δu` is then a path sum over `≥ 2` hinge edges of `Γ`
   (T1 Step 8) and `2Δu = (σ_g − σ_f)x_{v_w}` is a non-trivial affine condition on `X`; it can hold
   at isolated `t`, and nothing in H1–H6 forbids it holding identically. **[R2 correction: D-L1-3]**
   It is excluded only by H-L3's second clause. Measured: **0** occurrences.
7. **[N]** Check L1, block "L1.3": of `30 119 368` candidates, `239 892` have `h ≡ 0`; **all
   239 892** are coincident-copy pairs whose shared edge is a **hinge** with `src(e) = v_w`
   (type (a)); **0** are split-adjacent, **0** have no shared edge, and **0** have
   `v_w ∉ {v_a, v_b}` (the `λ ∉ {0,1}` branch of (L1.9) is empty on this corpus).
   Sharper still: **all 239 892** have `w == a` or `w == b` **as `M′`-vertices**, not merely as
   copies of the same `M`-vertex — the hinge point is welded by `cut.hpp::corner_to_prime` into a
   single `M′`-vertex shared by the two faces. So the identity class is decided by an **index
   comparison**, the cheapest possible combinatorial filter, and that is exactly what
   `contact.cpp`'s `if (p == a || p == b) continue;` already does, *before* the numeric identity
   test of T3.H.1. On this corpus the numeric identity test therefore fires on nothing the index
   comparison has not already removed. ∎

**What it replaces.** core.md's `[D]` paragraph in T5.2b.2 ("The excluded case is `A = C = B = 0` …
they are the **permanent incidences** of T3.H.1/T3.H.2 — a hinge point is a vertex of both incident
faces at every `θ`") and `check.md` R3.1a's "On the `squares` case I checked what they are: all 144
of them have `w` coincident with `a` or with `b` at every `θ` … 144/144, none of any other kind."
Both are assertions about one pattern; L1.3 proves the characterisation and Check L1 measures it on
84 graphs and 239 892 instances — **under H-L3 as a theorem, and unconditionally as a corpus
measurement [R2 correction: D-L1-3]**. It also settles the `λ`-interior case that the spec asks about
(`a vertex sitting on an edge it is incident to`): it is *possible* — condition 3 with
`λ ∈ (0,1)` — and it does not occur on this corpus, because a hinge point interior to another
face's edge requires three collinear flat vertices, which H4's simple positively-oriented faces do
not produce here.

## Corollary L1.4 — the `NOROOT` atom list is a combinatorial classification plus a measure-zero remainder

**Statement.** Fix a pattern `(M, σ)` and let `𝒞` be the ordered candidate list. Then:

**(i)** The partition of `𝒞` into `{h ≡ 0}` ⊔ `{class 3}` ⊔ `{class 2}` ⊔ `{class 1}` is determined
by `t` only through the two polynomials `h(0)` and `h′(0)`, and on the **coincident-copy** sublist
it is determined by the single vector `dS(t)` (Lemma L1.2). In particular:

* the `h ≡ 0` sublist is **fixed once per pattern**, by (L1.9) — no numeric test in `t` — **under
  H-L3**; without H-L3 the split branch (L1.10) adds the condition `Δu(t) = 0`, which is linear in
  `t` and therefore *not* combinatorial (Lemma L1.3, **[R2 correction: D-L1-3]**). On every design
  that passes the `0⁺` test H-L3 holds and the sublist is decided by the `M′`-vertex index
  comparison `w == a ∨ w == b`;
* the coincident-copy sublist is fixed once per pattern, by the four types of L1.1(ii);
* within it, class 3 versus class 2 is the sign of one quadratic polynomial `det(d, dS)(t)`;
* every remaining pair is class 1 **off a proper algebraic subset of the shape space** — under
  H-L1 — and the certificate needs the numeric test `|h(0)| ≤ tol` only for those.

**(ii)** Hence the atom list `(T5.1c)/(T5.1d)/(T5.1e′)` can be assembled by combinatorial
classification of pairs plus a numeric test confined to the accidental set, and the accidental set
has Lebesgue measure zero in `t`.

**(iii)** *Honest caveat.* (i)–(ii) are proved; H-L1, on which the "measure zero" of (ii) rests, is
**false on the measured corpus** for 0.579 % of candidate-harmonic *instances* (L1.1(iv); the unit
is instances, not distinct pairs — **[R2 correction: D-L1-4]**). For those pairs `h(0) ≡ 0`
structurally without being coincident-copy pairs, so a purely combinatorial classifier
under-counts the `C = 0` class. The consequence is *not* unsoundness: routing such a pair to
(T5.1c) when it belongs in (T5.1d) is exactly the failure core.md T5.2b.2's caveat (1) describes,
and it is the **unsafe** direction — (T5.1c)'s first atom `g(0)·g(T) > 0` fails and the pair is
reported as having a root, which makes `NOROOT` *conservative*, i.e. `R(ε)` shrinks. So the
practical rule stands: **classify by `|h(0)| ≤ tol` numerically, and use the combinatorics only to
know that the answer is structural.** That is what `contact.cpp` already does.

**Proof.** (i) is L1.1(i)+(ii), L1.2 step 5 and L1.3; (ii) is L1.1(iii); (iii) is L1.1(iv) plus the
direction of failure of (T5.1c) recorded in core.md T5.2b.2. ∎

## Check L1

**Program** `derivations/scratch/check_l1.cpp`. **Build line** (in the file header; the repo's
`clang++` needs `DEVELOPER_DIR` pointed at Xcode on this machine, see the note at the end of this
section):

```
clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
  -I code/src -I code/apps derivations/scratch/check_l1.cpp \
  code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
  code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
  code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
  -o /tmp/check_l1 && /tmp/check_l1 140 20
```

**Corpus.** The 8 Phase-2 reference cases (`kill_common.hpp::reference_cases()`) plus
`make_graph(id, 18, 46, 220)` for `id = 0 … 139` (the K1a population: Voronoi / Delaunay /
quad-random with `assign_orientation_relaxation` σ), of which **84 graphs** survive projection,
deployability (`1e−7`) and positive orientation. Per graph: `X₀` plus up to 19 Gaussian null-space
perturbations at `0.03 · scale`, all candidate pairs of all ordered face pairs (no broad phase).

**Output** (`/tmp/check_l1 140 20`):

```
check_l1 -- graphs used: 84  (random requested 140, samples/graph 20)
  candidate harmonics examined            : 30119368

L1.1  combinatorial predicate  C == 0  <=>  vw in {va, vb}
  coincident-copy pairs (vw in {va,vb})   : 801816
  numerically |C| <= tol                  : 999537
  AGREEMENTS                              : 29921647 / 30119368
  shared but C != 0  (must be 0)          : 0
  C == 0 but NOT shared (flat-collinear)  : 197721
     ... on rigid graphs (dim_null = 0, one sample, undecidable) : 22251
     ... of these, persistent over every sample of the graph : 174354
     ... of these, all three source vertices frozen          : 135410

L1.2/L1.3  class of a coincident-copy pair from dS alone
  numeric classes over ALL pairs: zero 239892  class1 29119831  class2 758176  class3 1469
  predicted == numeric on coincident-copy pairs : 801816 agree, 0 disagree
  closed forms  max |p - sf<d,dS>/2| = 1.990e-12   max |q + p| = 1.990e-12
                max |r - det(d,dS)/2| = 3.411e-13
  hinge dS = 2 sigma_g J (x_src - x_v):  479784 pairs, max err 1.000e-13
  class 3 occurrences: persistent 58   transient (isolated t) 152
                       undecidable (rigid graph, 1 sample) 1259

L1.3  the identically-zero pairs
  zero & coincident-copy                  : 239892
     ... of which w IS the same M'-vertex as a or b : 239892
     shared edge is a HINGE with src == v : 239892
     shared edge is a SPLIT (Delta u = 0) : 0
     no shared edge between f and g       : 0
  zero with vw NOT in {va,vb} (lambda interior) : 0

  coincident-copy pairs by shared-edge type: hinge 479784  split 159816  none 162216
  class-3 pairs by shared-edge type: hinge 160  split 24  none 100
  class-3 pairs across a HINGE: 160, of which <d, x_src - x_v> = 0 : 160
  split dS constant over the two endpoints (T1.B): 159816 pairs, max err 8.644e-14
  persistent non-shared C==0, by #frozen source vertices: 0:18762  1:20182  2:0  3:135410
```

Note on the class-3 line: `160 + 24 + 100 = 284` of the `1 469` class-3 occurrences are
coincident-copy pairs (the ones Lemma L1.2 covers). The other `1 185` are **not** coincident-copy
pairs — they are flat-collinear triples of three distinct vertices that additionally have `r = 0`,
i.e. they sit in the H-L1 remainder of L1.1(iv), and L1.2 says nothing about them beyond
(L1.2a), which applies to any class-3 harmonic whatever its combinatorial type.

**Reading of the numbers.**

| claim | measured |
|---|---|
| L1.1(ii) coincident copy `⟹ C = 0` | `801 816 / 801 816`, **0** counterexamples |
| L1.2 trichotomy from `dS` alone | `801 816` predictions, **0** disagreements |
| L1.2 closed forms `p, q, r` | `≤ 1.99e−12` (scale `~10²`) |
| L1.2b hinge `dS` formula | `479 784` pairs, `≤ 1.00e−13` |
| L1.2b hinge class 3 `⟺ d ⟂ e` | `160 / 160` |
| T1.B split `dS` endpoint-independent | `159 816` pairs, `≤ 8.64e−14` |
| L1.3 characterisation | `239 892 / 239 892` are hinge-`src` permanent incidences, and all of them have `w == a` or `w == b` as `M′`-vertices |
| H-L1 | **fails**: `174 354 / 30 119 368` structural non-coincident `C ≡ 0` |
| accidental `C = 0` | `1 116 / 30 119 368` = `3.7e−5` |

**Disagreements.** Zero, in every direction the lemmas actually claim. The one place where a
combinatorial rule and the numerics part company is H-L1, and that is stated as a **result**
(L1.1(iv)), not swept into a tolerance: the `197 721` cases are all of the form "three *distinct*
source vertices collinear in the flat state", `135 410` of them because all three vertices are
pinned by the fixed boundary and `38 944` because `tiling_3_4_3_12`'s shape space preserves a
collinearity among free vertices. Perturbing `t` does **not** remove them, which is the test the
spec asks for; only the residual `1 116` are removed by perturbation, and those are the accidental
set of L1.1(iii).

**Toolchain note (not a result).** On this machine `xcode-select` points at an Xcode whose
`xcodebuild` fails to load, so a bare `clang++` aborts. Both scratch programs were built with
`DEVELOPER_DIR=/Applications/Xcode.app/Contents/Developer` and the Xcode toolchain compiler
`.../Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++` with
`-isysroot .../MacOSX.sdk`. Compiler flags are otherwise exactly the header build line.

---

# Part L2 — the periodic dimension formula

## L2.0 Set-up: from the null vector to the `2 × 2` matrix

Throughout Part L2, `(M, σ)` is a **periodic** pattern given by a fundamental domain (one face per
translation class) together with the period matrix `P₀ := T = [t_h  t_v] ∈ R^{2×2}` (H-L2), and the
quotient system of `periodic_jacobian.hpp` supplies `X(t) = X₀ + Φ t` with `Φ ∈ R^{n_q × k}`,
`k = dim_null`. `Φ` is a null-space basis of the **scalar** system `[L; e_pin]`, and the design
coordinate is `t ∈ R^{2k}`: `Φ` acts on the `x` and the `y` coordinate **separately**, and the code
orders `t` as `(t_{2i}, t_{2i+1}) = ` (`φ_i` on `x`, `φ_i` on `y`) — `kill_k7.cpp::shape_point`.

**[A] L2.0a (the period law for the face potential).** For a lattice vector `τ ∈ {t_h, t_v}` there
is a **constant** `w_τ ∈ R²` with

```
        u_{f+τ}  =  u_f  +  w_τ  +  σ_f τ / 2        for every face f .               (L2.1)
```

*Proof.* Check that the right-hand side satisfies T1's increment law (T1.3) on every hinge edge of
the translated picture. Across `e + τ` between `f + τ` and `g + τ` the required increment is
`σ_g x_{src(e)+τ} = σ_g(x_{src(e)} + τ)`, while the ansatz gives
`(u_g + w_τ + σ_g τ/2) − (u_f + w_τ + σ_f τ/2) = σ_g x_{src(e)} + (σ_g − σ_f) τ/2 = σ_g x_{src(e)} + σ_g τ`,
using `σ_f = −σ_g` on a hinge (0.2). The two agree. Since `Γ` is connected (H2), a solution of the
increment system is unique up to one additive constant, so `w_τ` exists and is a single vector,
independent of `f`. ∎ This is the header derivation of `periodic_jacobian.hpp`, written out.

**[A] L2.0b (the deployed period).** From (T1.5) and (L2.1),

```
   y_{(v+τ, f+τ)}(θ) − y_{(v,f)}(θ)
        = c τ + s J ( 2u_{f+τ} − 2u_f − σ_f τ )  =  c τ  +  2 s J w_τ ,               (L2.2)
```

the `σ_f τ` terms cancelling. So with `Q := 2 J [w_h  w_v]`,

```
        P_θ = c P₀ + s Q ,     J(θ) = P_θ P₀⁻¹ = c I + s K ,     K = Q P₀⁻¹ ,         (L2.3)
```

`K` **constant in `θ`** — `[F, F33 / KILL_REPORT §K7 C1]`, measured `1.37e−13` against an
independent forward-kinematics call.

**[D] L2.0c (the period-potential increment of a null vector).** Fix a spanning tree of `Γ` and a
base face `f₀`. By (T1.6) and (L2.1),

```
        w_τ  =  Σ_{e ∈ path(f₀ → f₀+τ)} σ_{head(e)} x_{src(e)}  −  σ_{f₀} τ / 2 ,     (L2.4)
```

which is **affine** in `X` (linear part from the path sum, constant part `−σ_{f₀}τ/2`, since the
periods are fixed data). For a scalar function `φ ∈ R^{n_q}` on quotient vertices define

```
        d_τ(φ)  :=  Σ_{e ∈ path(f₀ → f₀+τ)} σ_{head(e)} φ_{src(e)}   ∈ R .            (L2.5)
```

`d_τ` is a linear functional on `R^{n_q}`; its restriction to `ker L` is independent of the path
(hence of the tree), because a change of path adds a cycle of `Γ`, and the closure of `δ` on cycles
is exactly `L φ = 0` (T1 Step 5/6). The matrix `D ∈ R^{k×2}` of `AchievableSet` is

```
        D_{i,1} = d_{t_h}(φ_i) ,     D_{i,2} = d_{t_v}(φ_i) ,     row i =: d_iᵀ .     (L2.6)
```

**[A] L2.0d (coordinate separation).** The linear part of (L2.4) is a signed sum of **position
vectors** with coefficients in `{0, ±1}` that do not depend on the coordinate. Hence for a
perturbation `δX = φ ⊗ e` with `e ∈ R²` fixed,

```
        δ w_τ  =  d_τ(φ) · e .                                                        (L2.7)
```

This is the whole reason the answer doubles: one scalar increment per null vector per period, and
the coordinate direction is carried along untouched.

## Lemma L2.1 — `dim 𝒦 = 2 · rank(D)`

**Statement.** Write `K(t) = K₀ + Σ_j t_j M_j` (affine, `[F, F33]`) and
`𝒦 = { K(t) : t ∈ R^{2k} }`, `A ∈ R^{4×2k}` the matrix whose columns are `vec(M_j)`, so
`dim 𝒦 = rank(A)`. Under H1–H3 and H-L2, with `e_x, e_y` the standard basis of `R²`:

```
   ┌──────────────────────────────────────────────────────────────────────────────┐
   │   M_{2i}   =   2 e_y d_iᵀ P₀⁻¹ ,        M_{2i+1}  =  − 2 e_x d_iᵀ P₀⁻¹ ,     │
   │   equivalently   M_{2i} P₀ = 2 e_y d_iᵀ  ,   M_{2i+1} P₀ = − 2 e_x d_iᵀ .    │
   │                                                                    (L2.8)    │
   │   dim 𝒦  =  rank(A)  =  2 · rank(D) .                             (L2.9)    │
   └──────────────────────────────────────────────────────────────────────────────┘
```

In particular `dim 𝒦 ∈ {0, 2, 4}` — never `1`, never `3`.

**Proof.**

1. **[A]** `t ↦ K(t)` is affine: `w_τ` is affine in `X` by (L2.4), `Q = 2J[w_h w_v]` is affine in
   `X`, `P₀` is *fixed* (the periods are lattice data, not design variables), and `X` is affine in
   `t`. So `K = Q P₀⁻¹` is affine, and `M_j = ∂K/∂t_j` is the constant matrix
   `2 J [δ_j w_h  δ_j w_v] P₀⁻¹` where `δ_j` is the derivative along the `j`-th design direction.
2. **[A]** *The `x`-copy.* `t_{2i}` moves `X` by `φ_i ⊗ e_x`. By (L2.7),
   `δ w_h = d_{t_h}(φ_i) e_x` and `δ w_v = d_{t_v}(φ_i) e_x`, so
   `[δw_h  δw_v] = e_x d_iᵀ` and `M_{2i} P₀ = 2 J e_x d_iᵀ = 2 e_y d_iᵀ`, since `J e_x = e_y`.
3. **[A]** *The `y`-copy.* `t_{2i+1}` moves `X` by `φ_i ⊗ e_y`, giving `[δw_h δw_v] = e_y d_iᵀ` and
   `M_{2i+1} P₀ = 2 J e_y d_iᵀ = − 2 e_x d_iᵀ`, since `J e_y = − e_x`. This proves (L2.8), and
   with it the two "vanishing row" statements: `M_{2i} P₀` has first row `0` and second row
   `2 d_iᵀ`; `M_{2i+1} P₀` has first row `− 2 d_iᵀ` and second row `0`.
4. **[A]** *The span.* A general element of `span{M_j}` is
   ```
        Σ_i t_{2i} M_{2i} + Σ_i t_{2i+1} M_{2i+1}
             = 2 ( e_y aᵀ − e_x bᵀ ) P₀⁻¹ ,     a := Σ_i t_{2i} d_i ,  b := Σ_i t_{2i+1} d_i .
   ```
   The even and odd coefficients are free and independent, so `(a, b)` ranges over
   `row(D) × row(D)` — the *whole* product, `row(D)` being the span of `{d_i}`.
5. **[A]** *Injectivity.* The map `Ψ(a,b) := 2( e_y aᵀ − e_x bᵀ ) P₀⁻¹` is linear, and
   `2(e_y aᵀ − e_x bᵀ)` is the `2 × 2` matrix with first row `−2 bᵀ` and second row `2 aᵀ`; it
   vanishes iff `a = b = 0`. Right multiplication by the invertible `P₀⁻¹` (H-L2) is a bijection of
   `R^{2×2}`. So `Ψ` is injective on `R² × R²`, a fortiori on `row(D) × row(D)`.
6. **[A]** Therefore `dim span{M_j} = dim(row(D) × row(D)) = 2 rank(D)`, and
   `rank(A) = dim span{vec(M_j)} = dim span{M_j}` because `vec` is a linear isomorphism
   `R^{2×2} → R⁴`. That is (L2.9). Since `D` has two columns, `rank(D) ∈ {0,1,2}`, so
   `dim 𝒦 ∈ {0,2,4}`. ∎

7. **[N]** Check L2: `133` patterns (the 33 K7 patterns plus 100 fresh random Voronoi tori),
   `rank(A) = 2 rank(D)` on **133/133**; the factorisation (L2.8) reconstructed from row `i` of `D`
   matches the numerically assembled `M_j` to `5.27e−14`, and the row of `M_j P₀` that (L2.8) says
   is zero is zero to `1.59e−14`. Observed `dim 𝒦 ∈ {0, 2, 4}` on all 133 — the value `2` occurs
   (`squares_3x2`), so the parity statement is not vacuous.

**What it replaces.** `results/kill/KILL_REPORT.md` §K7 C1's bullet "**`dim 𝒦 = 2·rank(D)` on
33/33**, where `D` is the `k × 2` matrix of period-potential increments `w_t` — this is the formula,
and it is exact", and `STATE.md` F33's "`dim𝒦 = 2·rank(D)` on 33/33". Both are measurements of a
formula that had no derivation. Lemma L2.1 derives it, and (L2.8) gives more than the dimension: the
**explicit generator**, so `K` can be written down from `D` without ever calling the forward
kinematics.

**Where `P₀` singular would break it.** Only step 5. If `P₀` were singular, `Ψ` could have a
kernel and `dim 𝒦` would drop below `2 rank(D)`; but `P₀ = [t_h t_v]` is the basis of a rank-2
lattice, so `det P₀ ≠ 0` always (H-L2). `[N]` `‖P₀ − T‖_∞ ≤ 6.22e−14` on 33/33
(`KILL_REPORT` §K7).

## Lemma L2.2 — what `rank(D)` is

**Statement.** Under H1–H3, H-L2:

**(a)** `rank(D) ≤ min(2, k)`, hence `dim 𝒦 ≤ 4`, and `dim 𝒦 ≤ 2k`. Both bounds are attained.

**(b)** `rank(D) = 0` **iff** every null vector has zero period-potential increment along **both**
periods, i.e. `d_{t_h}(φ) = d_{t_v}(φ) = 0` for every `φ ∈ ker L`; equivalently `K(t) ≡ K₀`, the
period Jacobian is **frozen** over the entire shape space, and `J(θ) = cos(θ/2) I + sin(θ/2) K₀`
is the same map for every design.

**(c)** *(the combinatorial meaning)* Reading `d_{t_h}, d_{t_v}` as covectors on all of `R^{n_q}`
via (L2.5),

```
        rank(D)  =  rank( [ L ; d_{t_h} ; d_{t_v} ] )  −  rank( L ) .                (L2.10)
```

So `rank(D)` counts how many of the two period-potential increments are **not already determined by
the deployability equations**: `rank(D) = 0` says both period covectors lie in the row space of `L`,
so Eq. (2) already fixes `w_h` and `w_v` and no design freedom can move them.

**(d)** *(the `squares_3x3` worked case, [N] + [A])* `squares_3x3` on the quotient has `k = 4` and
`D = 0` to `3.08e−18`, so `rank(D) = 0` and `dim 𝒦 = 0` although the shape space is
4-dimensional; `K₀` has `tr K₀ = 0.000000` and `det K₀ = 1.000000`, i.e. `K₀` is a rotation by
`π/2` and `J(θ) = cos(θ/2) I + sin(θ/2) K₀` is a **rigid rotation by `θ/2`** at every design and
every angle: the cell never opens, no matter what the four design coordinates do. By (c) the reason
is that both period covectors already sit in `row(L)`, which Check L2 verifies directly through
(L2.10). This is the pattern on which the project's earlier guess `min(4, 2·dim_null)` fails, and
Lemma L2.1 explains the failure precisely: **shape-space dimension does not bound designability;
the rank of the period-potential map does.**

**(e)** *(sufficient condition for `rank(D) = 2` — **CONJECTURE**, not proved)*
`rank(D) = min(2, k)` for every pattern except those whose symmetry forces `D = 0`. Evidence:
**132 of 133** patterns of Check L2 satisfy `rank(D) = min(2, k)` exactly; the sole exception is
`squares_3x3`. I could not find a proof, and I could not find a combinatorial predicate on the split
cuts that separates `squares_3x3` from `squares_3x2` (which has `k = 1` and `rank(D) = 1`); so this
is printed as a conjecture with its evidence, per the spec's instruction, and **not** as a lemma.

**Proof.**

1. **[A] (a)** `D ∈ R^{k×2}`, so `rank(D) ≤ 2` and `rank(D) ≤ k`. With (L2.9), `dim 𝒦 ≤ 4` and
   `dim 𝒦 ≤ 2k`. Attained: `hexagons_2x2` has `k = 4`, `rank(D) = 2`, `dim 𝒦 = 4`;
   `squares_3x2` has `k = 1`, `rank(D) = 1`, `dim 𝒦 = 2` (the `2k` bound is tight there).
2. **[A] (b)** `rank(D) = 0 ⟺ D = 0 ⟺ d_{t_h}(φ_i) = d_{t_v}(φ_i) = 0` for every basis vector,
   hence by linearity for every `φ ∈ ker L`. By (L2.8) all `M_j = 0`, so `K(t) ≡ K₀`. Conversely
   `K ≡ K₀` gives `A = 0` and, by (L2.9), `rank(D) = 0`.
3. **[A] (c)** Let `N := ker L ⊆ R^{n_q}`. The restriction map `(R^{n_q})* → N*` is linear and
   surjective with kernel the annihilator of `N`, which is `row(L)`. Hence
   `dim span{ d_{t_h}|_N , d_{t_v}|_N } = dim( span{d_{t_h}, d_{t_v}} + row(L) ) − dim row(L)`,
   which is the right-hand side of (L2.10). The left-hand side is `rank(D)`, since row `i` of `D`
   is `(d_{t_h}(φ_i), d_{t_v}(φ_i))` and the `φ_i` are a basis of `N`. ∎
4. **[N] (d)** Check L2, the `[rank(D) = 0 with k = 4]` line and the covector test.

**What it replaces.** `KILL_REPORT.md` §K7's sentence "`squares_3x3` has `dim_null = 4` (6 split
cuts) yet `dim 𝒦 = 0`: its `K` is frozen at a **rotation** … Shape-space dimension therefore does
**not** bound designability; the rank of the period-potential map does." That is a correct
observation with no mechanism; (b)+(c) supply the mechanism, and (e) says honestly that the
*positive* direction — predicting `rank(D) = 2` from the cut structure — is still open.

## Corollary L2.3 — conformal at `θ = 0` ⟹ conformal at every `θ`

**Statement.** `J(θ)` is a similarity (a conformal linear map) for **every** `θ` if and only if `K`
is one, which is the pair of **linear** equations in `t`

```
        K₁₁(t) = K₂₂(t) ,        K₁₂(t) = − K₂₁(t) .                                  (L2.11)
```

**Proof. [A]** The similarities of `R²` are `{ aI + bJ : (a,b) ≠ (0,0) }`, a 2-dimensional linear
subspace of `R^{2×2}` minus the origin, and (L2.11) is exactly the statement `K ∈ span{I, J}`. If
`K = αI + βJ` then `J(θ) = c I + s K = (c + sα) I + (sβ) J ∈ span{I,J}`, and it is non-singular
because `det J(θ) = (c + sα)² + (sβ)²`, which vanishes only if `c + sα = sβ = 0`; the map is then
the zero matrix, excluded since `det J(0) = 1` and `det J(θ)` is a non-negative quadratic in `(c,s)`
that vanishes at most at isolated `θ`. Conversely, if `J(θ₀) ∈ span{I,J}` for some `θ₀` with
`s(θ₀) ≠ 0`, then `s K = J(θ₀) − cI ∈ span{I,J}`, so `K ∈ span{I,J}`. Each `K_{ij}` is affine in
`t` by L2.1 step 1, so (L2.11) is two affine equations. ∎

**[N]** `[F, KILL_REPORT §K7 C2]` the two equations are consistent on all 25 patterns with
`dim 𝒦 ≥ 1` (residual `≤ 9.67e−17`) and at the solution the conformal distortion of `J(θ)` over
200 angles spanning `[0, π]` is `≤ 2.31e−14`.

**What it replaces.** The 2026 paper's §5.1 sentence "once the derivative of the conformal
distortion vanishes at `θ = 0`, the deployment remains conformal for all opening angles", asserted
*empirically*; and the corresponding one-line note in `periodic_jacobian.hpp`'s header, which states
the argument but is not a numbered result anywhere. With L2.1 it is a two-line consequence, and by
(L2.8) the two equations (L2.11) can be written directly in terms of `D` and `P₀`.

## Check L2

**Program** `derivations/scratch/check_l2.cpp`. Its population, quotient/super-patch pipeline and
`achievable()` are copied **verbatim** from `code/apps/kill_k7.cpp` (lines 36–243) so the numbers
compare directly with `results/kill/k7/k7_main.csv`; nothing in `code/` is modified.

**Build line** (in the file header; same `DEVELOPER_DIR` note as Check L1):

```
clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
  -I code/src -I code/apps derivations/scratch/check_l2.cpp \
  code/src/core/*.cpp code/src/method/*.cpp -o /tmp/check_l2 && /tmp/check_l2 100
```

**Corpus.** The 33 K7 patterns (7 tiling families × {2×2, 3×2, 3×3} = 21, plus the 12 Voronoi tori
`voronoi_torus_0_n20 … _11_n200`) **plus 100 fresh random Voronoi tori** of 12–45 sites at 100
distinct `σ` seeds (`quotient_sigma` max-cut on the quotient dual). 133 of 133 built and solved.

**Ranks** are dense `JacobiSVD` with threshold `1e-10 · σ_max` (floor `1e-13`).

**Output** (tail; the full per-pattern table is 133 rows):

```
squares_2x2                  k=  0  dimK=0 rank(A)=0 rank(D)=0  2rD=0 OK  fact_err=0.00e+00
squares_3x2                  k=  1  dimK=2 rank(A)=2 rank(D)=1  2rD=2 OK  fact_err=6.94e-16
    [rank(D) = 0 with k = 4]  max|D| = 3.084e-18   tr K0 = 0.000000  det K0 = 1.000000
squares_3x3                  k=  4  dimK=0 rank(A)=0 rank(D)=0  2rD=0 OK  fact_err=2.06e-18
hexagons_2x2                 k=  4  dimK=4 rank(A)=4 rank(D)=2  2rD=4 OK  fact_err=8.33e-16
...
voronoi_torus_11_n200        k=204  dimK=4 rank(A)=4 rank(D)=2  2rD=4 OK  fact_err=5.27e-14

check_l2 -- patterns evaluated: 133
  rank(A) == 2 rank(D):  133 / 133    failures: 0
  rank(D) distribution:   0 -> 8   1 -> 1   2 -> 124
  max |M_j - formula(D)|            : 5.271e-14
  max |vanishing row of M_j T|      : 1.592e-14
  L2.2 covector test  rank([L;delta_h;delta_v]) - rank(L) == rank(D):  94 ok, 0 mismatch
```

| claim | measured |
|---|---|
| L2.1 `rank(A) = 2 rank(D)` | **133 / 133**, 0 failures |
| L2.1 (L2.8) generator formula | `max |M_j − 2 e_y d_iᵀP₀⁻¹|` etc. `= 5.271e−14` |
| L2.1 (L2.8) vanishing rows of `M_j P₀` | `≤ 1.592e−14` |
| L2.1 `dim 𝒦 ∈ {0,2,4}` | 133 / 133 (values 0, 2 and 4 all realised) |
| L2.2(c) `rank(D) = rank([L;δ]) − rank(L)` | **94 / 94** patterns with `n_q ≤ 60` |
| L2.2(d) `squares_3x3` | `max|D| = 3.08e−18`, `tr K₀ = 0`, `det K₀ = 1` |
| L2.2(e) `rank(D) = min(2,k)` | **132 / 133** (exception: `squares_3x3`) |

**Distribution of `rank(D)`.** `0 → 8` patterns (seven of them trivially, `k = 0`: `squares_2x2`,
all three `triangles`, all three `kagome`; the eighth is `squares_3x3` with `k = 4`),
`1 → 1` (`squares_3x2`, `k = 1`), `2 → 124`. So `dim 𝒦 = 4` on 124 of 133 — consistent with K7's
"25 of 33" once the 100 extra Voronoi tori, which are all `rank(D) = 2`, are included.

**Failures.** None. The covector test (L2.10) was run only on the 94 patterns with `n_q ≤ 60`
because it costs one `K` evaluation per quotient vertex; on the larger tori it is not measured, and
that is stated as a limitation, not as a pass.

---

# For the Checker

Read this section first; it names the three steps I would attack, and the exact statements to
attack them with. Everything else in this file I believe is either `[A]` algebra you can redo on
paper in ten minutes or `[N]` a number you can reproduce with the two build lines.

## Weakest step 1 — H-L1 is false, and I claim the failure is *safe*

**The statement to attack.** Corollary L1.4(iii): "routing a structurally flat-collinear
non-coincident pair to (T5.1c) is the **conservative** direction, so `NOROOT` under-reports
validity rather than over-reporting it."

**Why it is weak.** L1.1(iv) measures `174 354` of `30 119 368` candidate pairs (0.58 %) whose
`h(0)` vanishes identically on the shape space **without** being coincident-copy pairs. My argument
that this is safe is one line: with `g(0) = 0` the first atom `g(0)·g(T) > 0` of (T5.1c) is false,
so `NOROOT` is false and the pair is reported as having a root. That argument is about the **exact**
predicate. The code decides class membership by `|C| ≤ 1e−11 · scale²`, and near that threshold the
pair can be routed either way. If a structural `C ≡ 0` pair lands just **above** the tolerance it is
treated as class 1 and (T5.1c) is applied with `g(0)` at round-off: `g(0)·g(T)` then has the sign of
noise, and a `NOROOT = true` verdict on a pair that does have a root at `θ = 0⁺` would be **unsound**
— the wrong direction. I did not measure the near-threshold population.

**How to attack it.** Take the `tiling_3_4_3_12` reference case (Check L1 prints eight explicit
witnesses, e.g. `f=0 g=9 va=0 vb=1 vw=24`). Perturb `t` so that `|h(0)|` lands within a decade of
`1e−11 · scale²`, and check whether `validity_certificate` at that `t` returns `noroot = true` while
the exact `contact_angles` list contains an angle in `(0, ε)`. If it does, L1.4(iii) is wrong and
the certificate needs `|h(0)| ≤ tol` to be a *hard* structural predicate, not a numeric one.

## Weakest step 2 — the orthogonal frame assumes `a` and `b` lie in the SAME face

**The statement to attack.** L1.0b: "`C_ab = d`, `S_ab = − σ_f J d`, so `{C_ab, S_ab}` is an
orthogonal basis of `R²`."

**Why it is weak.** Every one of L1.0c, L1.1, L1.2 and L1.3 is a computation in that frame; if the
frame is not orthogonal, `D₀ = det(C_ab,S_ab)` can vanish and the four-way split (L1.5) collapses.
The identity `S_ab = −σ_f J C_ab` is core.md T3.3 and it holds **only** when `a` and `b` are two
corners of a single face: the `2u_f` terms cancel in the difference precisely because both copies
carry the same `u_f`. My candidate enumeration takes `(a,b)` from `c.prime_faces[f]` as consecutive
corners, so it holds by construction *in my check*. It is a hypothesis about the certificate's
candidate list, not a theorem about it.

**How to attack it.** I read `code/src/method/contact.cpp`'s `validity_certificate` and its
`scan` lambda: it takes `a = E[i]`, `b = E[(i+1) % ne]` from `E = PF[fe] = prime_faces[fe]`, so
`a` and `b` **are** consecutive corners of one face today, and L1.0b applies to the current
candidate list. Re-check this whenever the list changes. If any candidate has `pa` and `pb` in different faces — for instance if
some future broad phase pairs an `M′`-vertex against a *shared segment* of two faces rather than a
single face's edge — then L1.0b fails on it and Lemmas L1.1(ii), L1.2 and L1.3 say nothing about
that pair. Also check the `T4.2` shared-segment branch of core.md T5.2b″ Case A, whose substitute
pair `π′` is constructed from the far-side edges of two *different* faces.

## Weakest step 3 — the coordinate separation (L2.7), and the non-circularity of the `D` check

**The statement to attack.** L2.0d: "for `δX = φ ⊗ e` with `e ∈ R²` fixed, `δ w_τ = d_τ(φ) · e`,
with the **same scalar** `d_τ(φ)` for `e = e_x` and for `e = e_y`."

**Why it is weak.** The whole factor of two in `dim 𝒦 = 2 rank(D)` is this sentence. It relies on
`u` being a signed sum of *position vectors* with coordinate-independent coefficients (T1.6) and on
`Φ` being the null space of a **scalar** operator applied to each coordinate. The periodic quotient
system is `[L; e_pin] X = [−(T n_k)ᵀ ; x_pin]`: the right-hand side genuinely *differs* between the
two coordinates. I claim that this affects `X₀` only and not `Φ`, so the two copies of `φ_i` are
both null directions and carry the same `d_τ(φ_i)`. If that is wrong, `D` as computed by
`kill_k7.cpp` — which reads it off the **even** generators only (`if (j % 2 == 0)`) — is not the
object in Lemma L2.1 and the lemma is a coincidence.

**How to attack it.** The non-circular half of my check is the comparison of the **odd** generators
`M_{2i+1}` against row `i` of `D`, which `kill_k7.cpp` never touches; Check L2 reports
`5.271e−14` over 133 patterns. Attack it by breaking the symmetry deliberately: modify a scratch
copy of the pipeline so that the `x` and `y` copies of `φ_i` are scaled differently, or so that a
pattern has anisotropic periods with `T` far from a multiple of the identity (`voronoi_torus`
uses `T = L·I`; the tiling families do not, so `t3_4_3_12` and `snub_square` already test this), and
confirm the odd-generator error stays at `1e−14`. A second, independent attack: recompute `D`
directly from (L2.5) — a signed sum of `φ_i` over hinge sources along a path crossing the period —
and compare with `AchievableSet::D`. I did **not** do that; my `D` is the one the code computes.

## Two smaller things I would also probe

* **L2.2(e) is a conjecture, not a lemma.** `rank(D) = min(2, k)` on 132 of 133 patterns, with
  `squares_3x3` the exception. I have no combinatorial predicate that separates it from
  `squares_3x2`. If you find one, it belongs in L2.2 as part (e) proper.
* **L1.3's `λ`-interior branch is empty on this corpus, not empty in general.** Condition 3 of
  L1.3 with `λ ∈ (0,1)` describes a hinge point lying in the *interior* of another face's edge.
  Check L1 measures **0** such pairs on 84 graphs. Constructing one (a face with three collinear
  consecutive vertices, which H4 permits — simplicity does not forbid collinearity) would give a
  permanent incidence that is **not** a coincident-copy pair, and would show that the
  `h ≡ 0` sublist is not exhausted by L1.1(ii) type (a). That is the cleanest single
  counterexample available against this file, and I could not build it.

## Reproduction

```
# Check L1  (84 graphs, 30 119 368 candidate harmonics, ~0.6 s)
clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
  -I code/src -I code/apps derivations/scratch/check_l1.cpp \
  code/src/core/mesh.cpp code/src/core/cut.cpp code/src/core/holes.cpp \
  code/src/core/kinematics.cpp code/src/core/orientation.cpp code/src/core/generators.cpp \
  code/src/core/tutte_auxetic.cpp code/src/core/collision.cpp code/src/core/optimize.cpp \
  -o /tmp/check_l1 && /tmp/check_l1 140 20

# Check L2  (133 periodic patterns, ~1.0 s)
clang++ -std=c++20 -arch arm64 -O2 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
  -I code/src -I code/apps derivations/scratch/check_l2.cpp \
  code/src/core/*.cpp code/src/method/*.cpp -o /tmp/check_l2 && /tmp/check_l2 100
```

On this machine `xcode-select` points at an Xcode whose `xcodebuild` aborts, so a bare `clang++`
fails before it reaches the source. Prefix both lines with
`DEVELOPER_DIR=/Applications/Xcode.app/Contents/Developer` and use
`/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++`
with `-isysroot .../Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk`. This is an environment
defect, not a property of the code; `code/tests` were not run because they go through the same
compiler and this file changes no code.
