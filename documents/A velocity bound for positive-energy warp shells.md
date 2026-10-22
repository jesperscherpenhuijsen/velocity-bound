# A velocity bound for positive-energy warp shells

Oct 20, 2026 · @Jesper Scherpenhuijsen

## Abstract

Fuchs et al. (2024) exhibited a stationary metric built from a static matter shell and a shift vector confined to the shell wall that satisfies all four pointwise energy conditions at a shift amplitude of 0.02 c, and left the upper limit on that amplitude to future work. We determine it. For any static shell and radial shift profile S(r), the Eulerian energy density and stresses change at second order in the amplitude v while the momentum density is linear in v, so the null condition reads ρ⁰ + S⁰(e,e) − 2v J·e + O(v²) and the threshold v\* = min (ρ⁰ + S⁰(e,e))/(2 J·e) is exact at first order. In the weak-field limit |J| = |S'' + S'/r|/16π on the equator; since a unit transition over width Δ has max|S''| ≥ 4/Δ², the threshold obeys v\* ≤ 2π(ρ⁰ + p\_t⁰)\_max Δ² (G = c = 1), and for a uniform shell v\* ≤ (3/4) C (1 − R₁/R₂)(1 + p\_t/ρ). The published shell's threshold is 0.031 c by bisection on the full energy conditions, against 0.0316 c from the first-order formula and 0.052 c from the bound; eight shells follow the (ρ + p\_t)Δ² scaling and stay under the bound. The eight shells computed reach 1% to 6% of c; the bound places relativistic amplitudes at compactnesses near the construction's own limits. The finalized manuscript is the LaTeX file velocity\_bound.tex with its PDF, delivered alongside this document; this doc is its working version.

## 1. Introduction

Warp drive spacetimes transport observers along geodesics between two points of an asymptotically flat background by modifying the metric in a bounded region that moves with them (Alcubierre 1994; Natário 2002). Every member of the Natário class, the metrics with unit lapse and flat spatial slices that contain the Alcubierre, Natário, Lentz and Fell and Heisenberg solutions, violates the null energy condition (Santiago, Schuster and Visser 2022). Bobrick and Martire (2021) argued that a physical warp drive must therefore carry a non-unit lapse or a curved spatial metric, positive ADM mass, and a subluminal velocity, and Fuchs et al. (2024) exhibited the first explicit example: a static, spherically symmetric matter shell of 4.49×10²⁷ kg between radii of 10 m and 20 m, with a shift vector added inside the shell wall, satisfying the null, weak, strong and dominant energy conditions at a shift amplitude of 0.02 c.

The same paper states the limitation that this work quantifies. The shift was kept below the level at which the momentum flux it induces causes energy-condition violations, and the authors note that an upper limit on the shift exists before momentum flux exceeds energy density, leaving its determination to future work. The question is therefore: how fast can a positive-energy warp shell move, and what sets the limit?

This paper answers both questions for the shell-plus-shift class.

1. For any static shell and any radial shift profile, the energy density and spatial stresses measured by Eulerian observers change only at second order in the shift amplitude v, while the momentum density is linear in v. The null energy condition along the null directions k = n + e therefore reads ρ⁰ + S⁰(e,e) − 2v J·e + O(v²), with J the linearized momentum density fixed by the momentum constraint. The threshold v\* = min over the shell of (ρ⁰ + S⁰(e,e))/(2 J·e) is exact at first order.
2. In the weak-field limit the linearized momentum density has a closed form in the shift profile S(r). At the equator |J| = |S'' + S'/r|/16π, so the threshold is set by the second derivative of the shift profile against the local energy density plus tangential pressure.
3. Because a monotone transition of unit height over a width Δ has max|S''| ≥ 4/Δ², the threshold obeys v\* ≤ 2π(ρ + p\_t)\_max Δ² up to a factor (1 − Δ/2R₁)⁻¹, in units with G = c = 1. For a uniform shell this becomes v\* ≤ (3/4) C (1 − R₁/R₂)(1 + p\_t/ρ), where C = 2GM/(R₂c²) is the compactness.
4. Numerically, the threshold of the Fuchs et al. shell is v\_NEC = 0.031 c by bisection on the full nonlinear energy conditions; the first-order formula gives 0.0316 c, the weak-field equatorial formula 0.036 c, and the bound 0.052 c. A scan over compactness, radius ratio and shift smoothing is reported in Section 5.

The practical content is a scaling law: the velocity of a positive-energy warp shell of this type is bounded by its compactness times its relative wall thickness. A shell with a wall half its radius needs C ≈ 0.27 to reach 0.1 c at the bound, and no shell with C below 2/3 reaches 0.5 c even when the bound is saturated, which places relativistic velocities in the regime where the shell itself is close to forming a horizon and where the weak-field expressions no longer apply.

## 2. Setup

Units are G = c = 1 unless SI values are written out; energy densities convert to SI by the factor c⁴/G. Signature (−,+,+,+). The metric is written in 3+1 form in the frame comoving with the shell, where it is stationary:

```latex
ds^2 = -\alpha^2 dt^2 + \gamma_{ij}\,(dx^i + \beta^i dt)(dx^j + \beta^j dt), \qquad g_{ti} = \beta_i = \gamma_{ij}\beta^j .
```

**The shell.** A static, spherically symmetric matter distribution with line element −e^{2a(r)}dt² + e^{2b(r)}dr² + r²dΩ². In Cartesian-type coordinates with n\_i = x\_i/r,

```latex
\alpha^0 = e^{a(r)}, \qquad \gamma_{ij} = \delta_{ij} + \left(e^{2b(r)} - 1\right) n_i n_j , \qquad \beta^0_i = 0 .
```

The shell carries a vacuum interior for r < R₁, matter for R₁ < r < R₂, and the Schwarzschild exterior with mass M and compactness C = 2M/R₂ for r > R₂. Its Eulerian stress-energy is diagonal: energy density ρ⁰(r), radial pressure p\_r⁰(r), tangential pressure p\_t⁰(r). Nothing below depends on how a and b are generated; Fuchs et al. build them from smoothed constant-density and Tolman-Oppenheimer-Volkoff profiles, and that construction is reproduced for the numerical part.

**The warp term.** Following equation (26) of Fuchs et al. (2024), a shift along the direction of motion x is added to the covariant time-space component and nothing else:

```latex
\beta_i = -\,v\, S(r)\, \delta_i^{\,x}, \qquad S(r) = 1 \ (r < R_1 + R_b), \quad S \to 0 \ (r > R_2 - R_b), \quad S \text{ monotone between},
```

with g\_tt left at −e^{2a}. The lapse therefore becomes α² = e^{2a} + v² S² γ^{xx}, which differs from α⁰ only at second order in v. The spatial metric does not change. Inside R₁ + R\_b the shift is constant, so the interior stays flat and the passengers are in vacuum; between R₁ + R\_b and R₂ − R\_b the shift changes, and that is where the energy conditions are tested. The analysis holds for any monotone S(r) and any static shell profile; the specific compact sigmoid and moving-average smoothing of Fuchs et al. enter only through the numbers in Section 5.

**Eulerian frame and energy conditions.** The Eulerian observer n^μ = (1/α, −β^i/α) and an orthonormal spatial triad e\_(i) built from γ\_ij give the frame components

```latex
\rho_E = T_{\mu\nu} n^\mu n^\nu, \qquad j_i = -\,T_{\mu\nu}\, n^\mu e_{(i)}^{\nu}, \qquad S_{ij} = T_{\mu\nu}\, e_{(i)}^\mu e_{(j)}^\nu ,
```

energy density, momentum density and spatial stress. Every null vector at a point is proportional to k = n + e for a unit spatial e, so the null energy condition is

```latex
T_{\mu\nu} k^\mu k^\nu = \rho_E - 2\, j_i e^i + S_{ij} e^i e^j \;\geq\; 0 \quad \text{for all unit } e .
```

The weak condition adds ρ\_E ≥ 0 and the same inequality for all timelike unit vectors; the dominant condition requires in addition that the energy flux be causal, |j| ≤ ρ\_E; the strong condition replaces T by T − ½T g in the weak inequality. For the stress tensors met here the null condition is the first to fail as v grows (Section 3), so it is the one tracked.

## 3. First-order expansion and the exact threshold

The metric depends on v only through β\_i = −vS δ\_i^x and through α² = e^{2a} + v²S²γ^{xx}. The 3+1 constraint and evolution equations of a stationary metric then sort the Eulerian stress-energy by its order in v.

**Energy density.** The Hamiltonian constraint,

```latex
16\pi\,\rho_E = {}^{(3)}\!R + K^2 - K_{ij}K^{ij},
```

contains the scalar curvature of γ\_ij, which does not depend on v, and the extrinsic curvature K\_ij = (D\_iβ\_j + D\_jβ\_i)/(2α), which is O(v). Hence ρ\_E = ρ⁰ + O(v²): the shift does not change the Eulerian energy density at first order.

**Momentum density.** The momentum constraint,

```latex
8\pi\, j_i = D_j\!\left(K^j{}_i - \delta^j_i K\right),
```

is linear in K\_ij, and K\_ij = v K̃\_ij /(1 + O(v²)) with

```latex
\tilde K_{ij} = \frac{D_i \tilde\beta_j + D_j \tilde\beta_i}{2\alpha^0}, \qquad \tilde\beta_i = -S(r)\,\delta_i^{\,x} .
```

Therefore j\_i = v J\_i + O(v³) with 8π J\_i = D\_j(K̃^j\_i − δ^j\_i K̃). The covariant derivative D is that of the shell's spatial metric and does not depend on v.

**Spatial stress.** For a stationary metric the evolution equation of K\_ij reads

```latex
0 = -D_iD_j\alpha + \alpha\left({}^{(3)}\!R_{ij} + K K_{ij} - 2K_{ik}K^k{}_j\right) + \mathcal{L}_\beta K_{ij} - 8\pi\alpha\left(S_{ij} - \tfrac{1}{2}\gamma_{ij}(S - \rho_E)\right).
```

Every term involving K or L\_β is O(v²), α and D\_iD\_jα are α⁰ and D\_iD\_jα⁰ up to O(v²), and ³R\_ij is v-independent. Hence S\_ij = S⁰\_ij + O(v²).

**The expansion.** Inserting these into the null condition of Section 2 gives, for every unit spatial e at every point x,

```latex
T_{\mu\nu}k^\mu k^\nu = \rho^0(x) + S^0_{ij}(x)\,e^i e^j - 2v\, J_i(x)\,e^i + O(v^2).
```

The zeroth-order part is non-negative wherever the static shell satisfies the null condition; the first-order part is negative for directions e with J·e > 0 and is unbounded below in v. Define at each point

```latex
v_1(x) = \min_{e:\ J\cdot e > 0}\ \frac{\rho^0(x) + S^0_{ij}(x)e^ie^j}{2\, J_i(x) e^i}, \qquad v^* = \min_{x \in \text{shell}} v_1(x).
```

**Theorem (first-order threshold).** To first order in the shift amplitude, the shell-plus-shift metric satisfies the null energy condition everywhere if and only if v ≤ v\*. The first violation appears at the minimizing point x\* and in the minimizing direction e\*, and the violated quantity grows linearly as 2(v − v\*)J·e\* for v slightly above v\*.

**Remarks.** (i) The second-order corrections are of relative size v, so for v\* of a few percent the first-order threshold is accurate to a few percent; Section 5 finds 2% for the Fuchs et al. shell. (ii) The dominant condition's flux inequality |j| ≤ ρ\_E fails at v\_DEC = min ρ⁰/|J| + O(v²). Since ρ⁰ ≥ S⁰(e,e) wherever the static shell obeys the dominant condition, v\_DEC ≥ v\*: the null condition is the binding one, and the weak and strong conditions, which contain the null one, fail together with it. (iii) The minimizing direction is close to the direction of J itself; replacing the minimum over e by e = Ĵ changes v\* by under 0.1% in the numerical cases below, because the tangential stress of the shell is nearly isotropic in the plane transverse to the radial direction. (iv) Nothing in the argument uses spherical symmetry of the shift; it uses only that the shift enters g\_ti linearly and leaves γ\_ij unchanged, and that the background is stationary.

## 4. Weak-field closed form and the compactness-thickness bound

In the limit of small compactness the shell's spatial metric is flat and its lapse is unity, so K̃\_ij = ½(∂\_iβ̃\_j + ∂\_jβ̃\_i) with β̃ = (−S(r), 0, 0). With n\_i = x\_i/r and a prime for d/dr,

```latex
\tilde K_{xx} = -S' n_x, \qquad \tilde K_{xy} = -\tfrac{1}{2} S' n_y, \qquad \tilde K_{xz} = -\tfrac{1}{2} S' n_z, \qquad \tilde K = -S' n_x ,
```

and the other components vanish. Inserting into 8πJ\_i = ∂\_j K̃\_ji − ∂\_i K̃ and using ∂\_i n\_j = (δ\_ij − n\_i n\_j)/r gives, with θ the angle from the direction of motion (cos θ = n\_x),

```latex
8\pi J_x = -\tfrac{1}{2}\left[ S'' \sin^2\theta + \frac{S'}{r}\left(1 + \cos^2\theta\right)\right], \qquad 8\pi J_y = \tfrac{1}{2}\, n_x n_y \left(S'' - \frac{S'}{r}\right), \qquad 8\pi J_z = \tfrac{1}{2}\, n_x n_z \left(S'' - \frac{S'}{r}\right).
```

Two locations are transparent. On the axis (θ = 0) the momentum density is radial with 8π|J| = |S'|/r, and the first-order threshold involves the radial pressure: v₁ = 4πr(ρ⁰ + p\_r⁰)/|S'|. On the equator (θ = π/2) the momentum density points along the direction of motion, which is tangential to the shell there, with

```latex
|J| = \frac{|S'' + S'/r|}{16\pi}, \qquad v_1^{\text{eq}}(r) = \frac{8\pi\left(\rho^0(r) + p_t^0(r)\right)}{\left| S''(r) + S'(r)/r \right|} .
```

For a transition of width Δ small compared with the radius, S'' is of order 1/Δ² while S'/r is of order 1/(Δ r), so the equatorial threshold is the smaller one and is governed by the second derivative of the shift profile against the local energy density plus tangential pressure. Section 5 confirms that the minimizing point lies on the equator at the inner edge of the transition, where |S''| is largest and the smoothed density has not yet reached its plateau.

**Lemma (second derivative of a unit transition).** Let S be C¹ on \[r\_a, r\_b\] with S(r\_a) = 1, S(r\_b) = 0, S'(r\_a) = S'(r\_b) = 0, and let Δ = r\_b − r\_a and m = max|S''|. Then |S'(r)| ≤ m·min(r − r\_a, r\_b − r) and 1 = ∫|S'| dr ≤ mΔ²/4, so m ≥ 4/Δ², with equality for the piecewise-quadratic profile whose S' is triangular. The same inequality gives |S'| ≤ mΔ/2.

**Bound.** Let r\_m be a point where |S''| = m. Since v\* is a minimum over the shell, v\* ≤ v₁^eq(r\_m). The denominator there is at least m − |S'(r\_m)|/r\_m ≥ m(1 − Δ/2R₁), hence, in the weak-field limit,

```latex
v^* \;\leq\; \frac{8\pi\left(\rho^0 + p_t^0\right)_{\max}}{m\,(1 - \Delta/2R_1)} \;\leq\; \frac{2\pi\left(\rho^0 + p_t^0\right)_{\max}\Delta^2}{1 - \Delta/2R_1}, \qquad G = c = 1 .
```

The factor (1 − Δ/2R₁)⁻¹ is needed only when S' and S'' have opposite signs at r\_m. In every case computed below the minimum of v₁ sits at the inner edge of the transition, where S' and S'' are both negative and the factor is absent; the simpler form v\* ≤ 2π(ρ⁰ + p\_t⁰)\_max Δ² held in all eight runs, with the ratio of the two sides between 1.4 and 3.4.

**Uniform shell.** For a shell of uniform density between R₁ and R₂, ρ⁰ = 3M/(4π(R₂³ − R₁³)) and Δ ≤ R₂ − R₁. Using R₂³ − R₁³ ≥ (R₂ − R₁)R₂², the simple form of the bound becomes

```latex
v^* \;\leq\; \frac{3}{4}\, C \left(1 - \frac{R_1}{R_2}\right)\left(1 + \frac{p_t}{\rho}\right), \qquad C = \frac{2GM}{R_2 c^2} .
```

The dominant energy condition of the static shell requires p\_t ≤ ρ, and for this construction it fails before C reaches 0.7 (Section 5). So positive-energy warp shells of this type are confined to v ≲ (3/2)C(1 − R₁/R₂) at the very most, and to a few tenths of that in practice because the density is not uniform where the transition is steepest and p\_t/ρ is well below one for shells that satisfy the dominant condition with margin. Reaching 0.1 c with a wall half the radius needs C of order 0.3 at the bound; 0.5 c needs C above 2/3, where the shell itself is near the dominant-condition limit and the weak-field expressions used here have corrections of order C.

**Curvature corrections.** The exact first-order J of Section 3 differs from the weak-field form by factors e^{−a}, e^{−b} and the Christoffel terms of γ\_ij. Section 5 measures the ratio of the exact to the weak-field equatorial |J| as 1.01 at C = 0.1, 1.07 at C = 1/3, 1.16 at C = 0.5 and 1.41 at C = 0.7: the curvature raises the momentum density and lowers the threshold, so the weak-field bound stays an upper bound at these compactnesses.

## 5. Numerical results

The published shell reaches 0.031 c; eight shells follow v\* ∝ (ρ + p\_t)\_max Δ² and stay under the bound.

**Method.** Warp Factory's shell generator was ported line by line to Python, the Einstein tensor computed with 4th-order finite differences on the z = 0 plane of a five-plane slab (h = 0.2 m, the published spacing), and the energy conditions evaluated in the Eulerian frame on 400 null directions and 10 rapidities plus a Hawking-Ellis eigenvalue test. Three quantities per shell: the first-order threshold v\* from J = j/v at v = 10⁻³ (linearity verified to 10⁻⁵); the weak-field equatorial formula from the radial profile; the true threshold by bisection on the full null condition, eight halvings, relative to the static shell's own numerical floor.

| Shell | C | R₁/R₂ | smoothing span | Δ \[m\] | ρ\_max \[J/m³\] | p\_t/ρ | v\* (first order) | v\_wf | v\_B / v\_U | v\_NEC (bisection) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| published | 0.333 | 0.50 | 4000 | 7.88 | 1.36×10⁴⁰ | 0.19 | 0.0316 | 0.0361 | 0.052 / 0.149 | 0.0305 to 0.0316 |
| lighter | 0.100 | 0.50 | 4000 | 7.88 | 4.07×10³⁹ | 0.04 | 0.0098 | 0.0101 | 0.014 / 0.039 | 0.0094 to 0.0105 |
| heavier | 0.500 | 0.50 | 4000 | 7.88 | 2.04×10⁴⁰ | 0.39 | 0.0456 | 0.0587 | 0.092 / 0.261 | 0.0422 to 0.0434 |
| heaviest (static DEC fails) | 0.700 | 0.50 | 4000 | 7.88 | 2.85×10⁴⁰ | 1.00 | 0.0601 | 0.0811 | 0.184 / 0.525 | 0.0527 to 0.0539 |
| thick wall | 0.333 | 0.25 | 4000 | 10.82 | 1.22×10⁴⁰ | 0.17 | 0.0529 | 0.0619 | 0.087 / 0.220 | 0.0492 to 0.0504 |
| thin wall | 0.333 | 0.75 | 4000 | 5.20 | 1.57×10⁴⁰ | 0.17 | 0.0170 | 0.0178 | 0.026 / 0.073 | 0.0164 to 0.0176 |
| sharper shift | 0.333 | 0.50 | 2000 | 6.94 | 1.38×10⁴⁰ | 0.34 | 0.0273 | 0.0317 | 0.046 / 0.168 | 0.0258 to 0.0270 |
| smoother shift | 0.333 | 0.50 | 8000 | 10.40 | 1.04×10⁴⁰ | 0.15 | 0.0404 | 0.0453 | 0.067 / 0.143 | 0.0387 to 0.0398 |

All velocities in units of c. v\_B = 2π(ρ + p\_t)\_max Δ²; v\_U = (3/4) C (1 − R₁/R₂)(1 + p\_t/ρ). The static shell satisfies the dominant condition in every row except C = 0.70, where p\_t reaches ρ.

**Reading the table.** The first-order threshold tracks the bisection value from above in every row, with a deficit growing from inside the bisection interval at C = 0.1 to 11% at C = 0.7, consistent with a second-order correction of relative size v\*. Relative to the published shell, (ρ\_max + p\_t,max)Δ² is 0.89, 1.28, 1.66 and 0.49 for the sharper shift, smoother shift, thick wall and thin wall, while v\* is 0.86, 1.28, 1.67 and 0.54 of its published value. Compactness at fixed radii scales v\* almost linearly (v\*/C from 0.098 to 0.086), the slow decline being the curvature correction to J (1.01 to 1.41). The bound holds in every row, by a factor 1.4 to 3.5 against the bisection midpoint; a mask 100 times lower changes neither v\* nor the bisection by more than 1%. The published shell's minimum sits on the equator at r = 12.4 m, 2.4 m into the wall, where S'' is most negative; max|S''| is 4.3% above the lemma's lower bound 4/Δ², so the published profile is nearly optimal already.

## 6. Discussion

The bound caps the shift amplitude of a positive-energy warp shell by the momentum density its wall can carry without violating the null condition, which is set by the wall's energy density. For the published shell that is 0.031 c; at ten percent compactness, 0.010 c. The profile factor is already within 4% of its optimum for the published smoothing, so the optimization Fuchs et al. anticipated can act only on (ρ + p\_t) where the shift is steepest and on Δ, both capped by compactness and wall thickness. A shift of 0.1 c needs C near 0.3 with a wall half the radius; 0.5 c with such a wall needs C above 2/3, where the construction's own dominant-condition margin has closed (it fails before C = 0.7); a thick-walled shell with p\_t = ρ would need C ≥ 1/3.

**Limitations.** First order in v (measured second-order corrections lower the threshold by up to 11%); weak-field closed form and bound with O(C) corrections measured, not derived (41% in |J| at C = 0.7); the class is the shell-plus-shift class with the shift confined to the wall. Whether the metrics correspond to a realizable matter model, and whether the shift amplitude is a transport velocity at all, are separate questions raised by Barzegar, Buchert and Vigneron (2026); the bound concerns the metric as written.

**Relation to the Natário-class no-go.** For unit lapse and flat slices every localized shift violates the null condition (Santiago, Schuster and Visser 2022). The shell adds a positive ρ⁰ and curved slices; Theorem 1 measures how far that escape lasts: while the first-order momentum term stays below ρ⁰ + p\_t⁰. As ρ⁰ → 0 the threshold goes to zero and the no-go is recovered.

**Physical reading.** A shift varying across the wall is a translational frame-dragging field, sourced by momentum density of order v/16πΔ². The null condition caps the momentum density a given energy density can carry at (ρ + p\_t)/2 along the flux direction. The threshold is the amplitude at which the wall would have to carry momentum comparable to its own energy density.

## 7. Conclusion

Adding a radial shift to a static shell changes the Eulerian momentum density at first order and everything else at second order, so the null condition fails at v\* = min (ρ⁰ + S⁰(e,e))/(2 J·e), with J set by the second derivative of the shift profile. Since that derivative is at least 4/Δ², v\* ≤ 2π(ρ + p\_t)\_max Δ², which for a uniform shell is (3/4) C (1 − R₁/R₂)(1 + p\_t/ρ). The published positive-energy warp shell reaches 0.031 c; eight shells follow the scaling and stay under the bound. The shells computed here reach 1% to 6% of c, and relativistic amplitudes require compactnesses near the construction's own dominant-condition and horizon limits.

## Appendix. Numerical method and validation

The pipeline is the one validated in the Phase 0 document: Einstein tensor within 1.5% of an exact symbolic Alcubierre result, Eulerian density within 7×10⁻⁶ of the closed form, Schwarzschild vacuum residual 4×10⁻⁴, Minkowski exactly zero. For the shell, halving the grid spacing at v = 0.2 c leaves the null minimum, its location and the frame tensor unchanged to three digits, and the static shell reproduces the 10³⁴ J/m³ numerical floor reported by Fuchs et al. The generator reproduces their steps: constant-density profile, constant-density TOV pressure, four moving-average passes (span 1.79 × sF on density, sF on pressure), mass function from the smoothed density, e^{2b} = (1 − 2Gm/rc²)⁻¹, a integrated inward from the Schwarzschild value, compact sigmoid shift smoothed twice, Cartesian γ\_ij = δ\_ij + (e^{2b} − 1) n\_i n\_j with g\_tx = −vS. Scripts: vbound.py (thresholds, bound, bisection), wfcheck.py (pipeline), vbound\_results.json (raw numbers), all delivered with the LaTeX manuscript.

## References

- Alcubierre M. 1994, The warp drive: hyper-fast travel within general relativity, Class. Quantum Grav. 11, L73.
- Barzegar H., Buchert T., Vigneron Q. 2026, General formalism, classification, and demystification of the current warp-drive spacetimes, [arXiv:2602.16495](https://arxiv.org/abs/2602.16495).
- Baumgarte T.W., Shapiro S.L. 2010, Numerical Relativity, Cambridge University Press.
- Bobrick A., Martire G. 2021, Introducing physical warp drives, Class. Quantum Grav. 38, 105009.
- Fell S.D.B., Heisenberg L. 2021, Positive energy warp drive from hidden geometric structures, Class. Quantum Grav. 38, 155020.
- Fuchs J., Helmerich C., Bobrick A., Sellers L., Melcher B., Martire G. 2024, Constant velocity physical warp drive solution, Class. Quantum Grav. 41, 095013, [arXiv:2405.02709](https://arxiv.org/abs/2405.02709).
- Helmerich C., Fuchs J., Bobrick A., Sellers L., Melcher B., Martire G. 2024, Analyzing warp drive spacetimes with Warp Factory, Class. Quantum Grav. 41, 095009; [code](https://github.com/NerdsWithAttitudes/WarpFactory).
- Lentz E.W. 2021, Breaking the warp barrier: hyper-fast solitons in Einstein-Maxwell-plasma theory, Class. Quantum Grav. 38, 075015.
- Natário J. 2002, Warp drive with zero expansion, Class. Quantum Grav. 19, 1157.
- Santiago J., Schuster S., Visser M. 2022, Generic warp drives violate the null energy condition, Phys. Rev. D 105, 064038.
