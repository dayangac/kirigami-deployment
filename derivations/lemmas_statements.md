# derivations/lemmas_statements.md — STATEMENTS ONLY (for the blind Checker pass)

Extracted mechanically by the orchestrator from derivations/lemmas.md: the standing hypotheses, the [D] definition paragraphs of the set-up sections L1.0 and L2.0 (their [A] algebra is withheld), and each lemma's **Statement** block. Proofs, 'What it replaces', checks and the Deriver's weakest-step notes are omitted on purpose.

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
   (d) f, g share the vertex v_w but no edge             copies around a cut cone
```

**(iii)** *(converse, under H-L1)* If `v_a, v_b, v_w` are three distinct vertices, then either
`h_π(0)` vanishes identically on `𝕏` — in which case the three flat points are collinear for
**every** design — or `{ t : h_π(0; t) = 0 }` is the zero set of a non-zero polynomial of degree
`≤ 2`, hence a proper algebraic subset of `R^{2k}` of Lebesgue measure zero: an **accidental**
degeneracy. H-L1 is exactly the assumption that the first alternative never occurs.

**(iv)** *(H-L1 fails, and how)* On the corpus of Check L1 (84 graphs, 30 119 368 candidates)
**197 721** candidates (0.66 %) have `|h_π(0)| ≤ tol` with `v_w ∉ {v_a, v_b}`. Of those,
**174 354** vanish at **every** shape-space sample of their graph, so they are structural, not
accidental: **135 410** have all three source vertices *frozen* (`Φ`-row zero — collinear triples
pinned by the fixed boundary), and **38 944** have one or two free vertices that the shape space
nevertheless keeps collinear (every witness the program printed is in `tiling_3_4_3_12`; I did not
enumerate the rest, so "all in `tiling_3_4_3_12`" is not claimed). A further **22 251** sit on graphs with
`dim_null = 0`, where there is no shape space to test. Only **1 116** of 30 119 368 (`3.7e−5`)
vanish at some samples and not others — the accidental set of (iii).

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

if `f` and `g` are **split-adjacent**, it reads `Δu = 0`, i.e. the two faces carry the same face
potential and the split cut never opens. So under H5 and the standing hypotheses the
identically-zero pairs are exactly

> the pairs in which `w` is the copy in `g` of the **hinge point** `src(e)` of an edge `e` shared
> by `f` and `g`, and that hinge point lies on the line of the tested edge `(a,b)` of `f` — in
> particular every pair of type (a) of L1.1(ii), where the hinge point *is* an endpoint of `(a,b)`.

This is a **combinatorial** description: it is decided by `(M, σ)` and the flat incidences, not by
`t`, except for the single collinearity in (L1.9), which for `λ ∈ {0,1}` is automatic.

## Corollary L1.4 — the `NOROOT` atom list is a combinatorial classification plus a measure-zero remainder
**Statement.** Fix a pattern `(M, σ)` and let `𝒞` be the ordered candidate list. Then:

**(i)** The partition of `𝒞` into `{h ≡ 0}` ⊔ `{class 3}` ⊔ `{class 2}` ⊔ `{class 1}` is determined
by `t` only through the two polynomials `h(0)` and `h′(0)`, and on the **coincident-copy** sublist
it is determined by the single vector `dS(t)` (Lemma L1.2). In particular:

* the `h ≡ 0` sublist is **fixed once per pattern**, by (L1.9) — no numeric test in `t` (L1.3);
* the coincident-copy sublist is fixed once per pattern, by the four types of L1.1(ii);
* within it, class 3 versus class 2 is the sign of one quadratic polynomial `det(d, dS)(t)`;
* every remaining pair is class 1 **off a proper algebraic subset of the shape space** — under
  H-L1 — and the certificate needs the numeric test `|h(0)| ≤ tol` only for those.

**(ii)** Hence the atom list `(T5.1c)/(T5.1d)/(T5.1e′)` can be assembled by combinatorial
classification of pairs plus a numeric test confined to the accidental set, and the accidental set
has Lebesgue measure zero in `t`.

**(iii)** *Honest caveat.* (i)–(ii) are proved; H-L1, on which the "measure zero" of (ii) rests, is
**false on the measured corpus** for 0.58 % of pairs (L1.1(iv)). For those pairs `h(0) ≡ 0`
structurally without being coincident-copy pairs, so a purely combinatorial classifier
under-counts the `C = 0` class. The consequence is *not* unsoundness: routing such a pair to
(T5.1c) when it belongs in (T5.1d) is exactly the failure core.md T5.2b.2's caveat (1) describes,
and it is the **unsafe** direction — (T5.1c)'s first atom `g(0)·g(T) > 0` fails and the pair is
reported as having a root, which makes `NOROOT` *conservative*, i.e. `R(ε)` shrinks. So the
practical rule stands: **classify by `|h(0)| ≤ tol` numerically, and use the combinatorics only to
know that the answer is structural.** That is what `contact.cpp` already does.

## L2.0 Set-up: from the null vector to the `2 × 2` matrix

Throughout Part L2, `(M, σ)` is a **periodic** pattern given by a fundamental domain (one face per
translation class) together with the period matrix `P₀ := T = [t_h  t_v] ∈ R^{2×2}` (H-L2), and the
quotient system of `periodic_jacobian.hpp` supplies `X(t) = X₀ + Φ t` with `Φ ∈ R^{n_q × k}`,
`k = dim_null`. `Φ` is a null-space basis of the **scalar** system `[L; e_pin]`, and the design
coordinate is `t ∈ R^{2k}`: `Φ` acts on the `x` and the `y` coordinate **separately**, and the code
orders `t` as `(t_{2i}, t_{2i+1}) = ` (`φ_i` on `x`, `φ_i` on `y`) — `kill_k7.cpp::shape_point`.

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

## Corollary L2.3 — conformal at `θ = 0` ⟹ conformal at every `θ`
**Statement.** `J(θ)` is a similarity (a conformal linear map) for **every** `θ` if and only if `K`
is one, which is the pair of **linear** equations in `t`

```
        K₁₁(t) = K₂₂(t) ,        K₁₂(t) = − K₂₁(t) .                                  (L2.11)
```


