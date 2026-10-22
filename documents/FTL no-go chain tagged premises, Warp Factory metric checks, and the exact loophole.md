# FTL no-go chain: tagged premises, Warp Factory metric checks, and the exact loophole

Oct 20, 2026 · @Jesper Scherpenhuijsen

## 1. Scope, conventions and tags

The chain below has one unproven link that matters: the achronal averaged null energy condition on self-consistent semiclassical backgrounds. Every other link is a theorem, a flat-space theorem, or a number this document recomputes.

**Conventions.** Signature (-,+,+,+); G = c = ħ = 1 unless SI units are written out. A warp metric is written in 3+1 form with lapse α, shift β^i and spatial metric γ\_ij. The Eulerian observer is the unit normal to the slices, n^μ = (1/α, -β^i/α); its energy density is ρ\_E = T\_μν n^μ n^ν. Energy conditions: NEC, T\_μν k^μ k^ν ≥ 0 for every null k; WEC, T\_μν u^μ u^ν ≥ 0 for every timelike u; DEC, WEC plus -T^μ\_ν u^ν causal; SEC, (T\_μν - ½ T g\_μν) u^μ u^ν ≥ 0. ANEC: ∫ T\_μν k^μ k^ν dλ ≥ 0 along a complete null geodesic with affine parameter λ. Achronal ANEC: the same, restricted to complete achronal null geodesics (no two points on the geodesic joined by a timelike curve).

**What counts as FTL.** Two definitions are used in the literature and both appear below. Olum (1998): a causal curve from region A to region B that arrives before every nearby null geodesic, defined on a flat neighbourhood at each end. Gao and Wald (2000): a null geodesic from past to future null infinity that arrives earlier than it would in the unperturbed spacetime, relative to the asymptotic structure. A spacetime passes the FTL test only if one of these holds; a large shift vector by itself does not.

**Tags.** Each premise carries exactly one tag.

| Tag | Meaning |
| --- | --- |
| PG | Proven in general: a theorem holding on arbitrary (or arbitrary globally hyperbolic) Lorentzian spacetimes, or a measurement whose interpretation does not depend on a background assumption |
| PF | Proven in flat spacetime: a theorem of quantum field theory or kinematics on Minkowski space, or a laboratory result whose analysis assumes a flat background |
| VN | Verified numerically: a number or inequality established by computation on specific metrics, in the cited paper or in Section 3 of this document |
| CJ | Conjecture: stated in the literature, with partial results, with no general proof in either direction |

A link is only as strong as its weakest premise. The chain's overall tag is the weakest tag among the premises that an FTL spacetime would have to defeat.

## 2. The theorem chain, link by link

Twelve links. Each gives the statement, the premises with tags, and what an FTL spacetime would have to defeat. Links 1 to 2 are kinematic, 3 to 7 are classical general relativity, 8 is quantum field theory, 9 to 12 are auxiliary.

### L1. Local speed limit

**Statement.** In any local inertial frame, no massive body and no signal carried by a known field moves faster than c.

**Premises.** P1.1 Local Lorentz invariance holds: cavity tests bound the anisotropy of c below 10^-17 (Nagel et al. 2015), gamma-ray burst timing pushes any linear energy dependence of the photon speed beyond the Planck scale (Abdo et al. 2009), and neutrino time of flight over 730 km is consistent with c (ICARUS 2012). Tag PF (experiment, flat-space analysis). P1.2 General relativity reduces to special relativity in every local frame (equivalence principle). Tag PG.

**What FTL must defeat.** Nothing here forbids a spacetime shortcut; L1 only removes the option of pushing a body past c through the local light cone.

### L2. Kinematic causality

**Statement.** If a device signals at u > c in its rest frame and obeys the relativity principle, two devices in relative motion v with u·v > c² return a reply before the question is sent (Tolman 1917; Benford, Book and Newcomb 1970). In general relativity the same closed loop appears for two warp bubbles in relative motion (Everett 1996) and for a wormhole whose mouths carry a relative time shift (Morris, Thorne and Yurtsever 1988).

**Premises.** P2.1 Lorentz transformations between inertial frames. Tag PF (experiment). P2.2 The device works identically in every inertial frame. Tag PG as a theorem premise; a warp bubble or wormhole evades it locally by singling out a frame, which is exactly why the GR versions of the loop need two devices. P2.3 Chronology protection: no closed timelike curve can form from regular initial data (Hawking 1992). Tag CJ. Partial support: for free fields in Hadamard states the renormalized stress tensor cannot stay regular everywhere on a compactly generated Cauchy horizon (Kay, Radzikowski and Wald 1997), tag PG within that setting.

**What FTL must defeat.** A working FTL device makes P2.3 load-bearing: either closed timelike curves are physical, or chronology protection destroys the second device. Neither branch is proven.

### L3. Superluminal metrics exist as solutions

**Statement.** The Alcubierre (1994), Natário (2002) and Krasnikov (1998) metrics and the Morris and Thorne (1988) wormhole are exact Lorentzian metrics. Each solves the Einstein equations for the stress tensor T = G/8π read off from its geometry.

**Premises.** None beyond the Einstein equations. Tag PG.

**What this does not say.** Existence of a metric says nothing about whether its T is realizable. Every later link is about T.

### L4. Olum's theorem

**Statement.** If a causal path between two flat regions is superluminal in Olum's sense, the null energy condition fails along the tangent of the fastest null geodesic, so the weak energy condition fails too (Olum 1998).

**Premises.** P4.1 Einstein equations. Tag PG. P4.2 Smooth Lorentzian metric with flat regions at both ends so that "faster than nearby light" is well defined. Tag PG (a definition, not a physical assumption).

**Tag of the link.** PG. The proof is the Raychaudhuri equation: a fastest null geodesic has no conjugate points, which forces R\_μν k^μ k^ν < 0 somewhere on it.

### L5. Perturbative superluminal censorship

**Statement.** In linearized gravity about Minkowski space, light cones can only be widened outward where the null energy condition fails (Visser, Bassett and Liberati 2000).

**Premises.** P5.1 Linearized Einstein equations about flat space. Tag PG within perturbation theory. P5.2 Asymptotic flatness. Tag PG.

**Tag of the link.** PG for weak fields. A strong-field FTL metric escapes L5 only to meet L4 and L6.

### L6. Gao and Wald: gravitational time delay

**Statement.** In a globally hyperbolic, asymptotically flat (or asymptotically AdS) spacetime in which the averaged null energy condition holds along every complete null geodesic and the null generic condition holds, no null geodesic from past to future null infinity arrives earlier than it would in the unperturbed spacetime (Gao and Wald 2000). Friedman, Schleich and Witt (1993) prove the companion result for topology: under ANEC every causal curve from past to future infinity is deformable to one that stays near infinity, so no wormhole shortcut exists.

**Premises.** P6.1 Global hyperbolicity. Tag CJ: it is the content of strong cosmic censorship, unproven in general. P6.2 Null generic condition. Tag PG in the sense that it fails only on measure-zero symmetric geometries such as exact plane waves. P6.3 ANEC along every complete achronal null geodesic of the spacetime, with T the expectation value of the quantum stress tensor in the actual state. Tag PF: proven in Minkowski space for unitary Lorentz-invariant quantum field theories (Faulkner, Leigh, Parrikar and Wang 2016; Hartman, Kundu and Tajdini 2017). On curved self-consistent backgrounds it is the Graham and Olum (2007) conjecture, tag CJ. Known partial results: proven for free fields perturbatively about flat space on classical backgrounds obeying the NEC (Kontou and Olum 2015); derived from the generalized second law, itself unproven (Wall 2010); violated along achronal geodesics on a conformally flat background that is not a self-consistent solution (Urban and Olum 2010); violated along non-achronal geodesics between Casimir plates.

**Tag of the link.** CJ, inherited from P6.3. This is the load-bearing link and the only one whose weakness an FTL proposal can exploit. Section 4 states it as a precise premise.

### L7. Generic warp drives violate the NEC

**Statement.** For every metric of the form ds² = -dt² + δ\_ij (dx^i - β^i dt)(dx^j - β^j dt) with a shift β that decays at infinity, the Eulerian energy density integrates to a non-positive total, is strictly negative somewhere unless β is a gradient, and the NEC fails for any localized non-trivial shift (Santiago, Schuster and Visser 2022). Appendix A rederives the integral statement in four lines.

**Premises.** P7.1 Unit lapse and flat spatial slices (the Natário class, which contains Alcubierre, Natário, Lentz and Fell and Heisenberg). Tag PG within the class. P7.2 Einstein equations. Tag PG.

**Tag of the link.** PG within the class. It settles the Lentz (2021) and Fell and Heisenberg (2021) positive-energy claims. It does not cover non-unit-lapse constructions (Bobrick and Martire 2021), which Section 3 checks numerically. Update: a lapse-weighted identity, ∫N(ρ\_E + S/3) = 2∫Nρ\_E + (1/12π)∮∇N·dS, now proves NEC violation for every stationary flat-slice drive with a one-component shift and any positive lapse, which covers the modified-time class (Scherpenhuijsen, draft 2026); the multi-component case remains a conjecture, tag CJ.

### L8. Quantum limits on negative energy

**Statement.** Quantum field theory allows negative energy density, but bounds its magnitude times duration (quantum inequalities) and forbids net negative energy along null lines in flat space (ANEC).

**Premises.** P8.1 ANEC in Minkowski space. Tag PF (as in P6.3). P8.2 Quantum null energy condition. Tag PF: free fields (Bousso, Fisher, Koeller, Leichenauer and Wall 2016), general interacting theories (Balakrishnan, Faulkner, Khandker and Wang 2019). P8.3 Quantum inequalities for free fields: flat space (Ford and Roman 1995, 1997; Fewster and Eveson 1998), tag PF; globally hyperbolic curved spacetimes with Hadamard states (Fewster 2000), tag PG for free fields. For interacting fields in four dimensions no quantum inequality is proven; two-dimensional conformal field theory is the only interacting case (Fewster and Hollands 2005). Tag CJ for interacting four-dimensional matter. P8.4 Application to warp metrics: with the flat-space inequality sampled on times short compared with the curvature radius, the Alcubierre wall is at most about 100 Planck lengths thick for v = c and the total negative energy for a 100 m bubble exceeds the mass-energy of the visible universe by about ten orders of magnitude (Pfenning and Ford 1997). Van Den Broeck (1999) shrinks the exterior and brings the figure to a few negative solar masses. The Krasnikov tube fares the same (Everett and Roman 1997), and the requirement stays unphysical at v much less than c (Lobo and Visser 2004). Tag VN for all four.

**Tag of the link.** PF for the inequalities that carry the numbers; VN for the numbers; CJ only for the extension to interacting matter, which can only tighten the bounds if the pattern of P8.1 and P8.2 holds.

### L9. Semiclassical pathologies of superluminal bubbles

**Statement.** A bubble with v > c carries a horizon pair. In the 1+1 reduction the renormalized stress tensor of a quantum field diverges at the front horizon as the bubble forms (Hiscock 1997; Finazzi, Liberati and Barceló 2009). Particles swept up at the front arrive inside blueshifted without bound as v → c (McMonigal, Lewis and O'Byrne 2012).

**Premises.** P9.1 Semiclassical gravity on the given background; 1+1 reduction for the divergence result. Tag VN (analytic in reduced dimension, not a four-dimensional theorem).

**Role.** Not needed for the no-go. It adds a stability requirement that any candidate must survive.

### L10. What has been made in the laboratory

**Statement.** Negative energy density exists relative to the vacuum between Casimir plates, with energy density -π²ħc/(720 d⁴), about -4×10^-4 J/m³ at d = 1 μm; the force has been measured (Lamoreaux 1997; Mohideen and Roy 1998). Squeezed vacuum states, measured to 15 dB (Vahlbruch et al. 2016), carry transient negative energy density bounded by P8.3. No direct laboratory measurement of a time-integrated negative energy density against a quantum inequality bound is known to me as of this writing. No negative gravitational mass has been observed.

**Premises.** P10.1 Casimir force measurements. Tag PF (experiment). P10.2 The inferred energy density follows from flat-space QFT. Tag PF.

**Role.** Fixes the gap between available and required negative energy at more than 40 orders of magnitude.

### L11. Coupling magnitude

**Statement.** G/c⁴ = 8.26×10^-45 m/J. One joule over one centimetre perturbs the metric by a fractional 10^-42; optical clocks resolve about 10^-20.

**Premises.** Measured constants. Tag PG.

**Role.** Independent of exotic matter: no laboratory energy produces a detectable metric change.

### L12. The wormhole route

**Statement.** A traversable wormhole throat requires NEC violation at or near the throat (Morris and Thorne 1988; Hochberg and Visser 1997, 1998). Creating one from trivial topology requires topology change, which in a compact region forces closed timelike curves or a singular metric (Geroch 1967) and violates the weak energy condition (Tipler 1977). Topological censorship (L6) then blocks the shortcut under ANEC.

**Premises.** P12.1 Flare-out geometry at the throat. Tag PG. P12.2 Lorentzian metric throughout the topology change. Tag PG. The ANEC premise is P6.3 again.

**Tag of the link.** PG for the throat and topology results; the shortcut itself inherits CJ from P6.3.

### Premise table

| Premise | Content | Tag | Where it is proven | Status if false |
| --- | --- | --- | --- | --- |
| P1.1 | Local Lorentz invariance | PF | Nagel et al. 2015; Abdo et al. 2009; ICARUS 2012 | L1, L2 fall |
| P2.3 | Chronology protection | CJ | Hawking 1992; Kay, Radzikowski and Wald 1997 | FTL implies closed timelike curves |
| P4.1, P7.2 | Einstein equations | PG | Classical GR | Whole chain replaced |
| P5.1 | Linearized gravity | PG (weak field) | Visser, Bassett and Liberati 2000 | Strong-field case falls to L4, L6 |
| P6.1 | Global hyperbolicity | CJ | Strong cosmic censorship, unproven | L6 inapplicable; L4 still holds |
| P6.2 | Null generic condition | PG (generic) | Gao and Wald 2000 | Only exact symmetric geometries escape |
| P6.3 | Achronal ANEC, self-consistent curved background | PF flat, CJ curved | Faulkner et al. 2016; Hartman et al. 2017; Graham and Olum 2007 | **The loophole** |
| P7.1 | Natário class (unit lapse, flat slices) | PG (in class) | Santiago, Schuster and Visser 2022 | Non-unit lapse needs separate check (Section 3) |
| P8.2 | QNEC | PF | Balakrishnan et al. 2019 | Local negative energy less constrained |
| P8.3 | Quantum inequalities | PF free fields; PG free curved; CJ interacting | Ford and Roman; Fewster 2000 | Interacting matter could evade magnitude bounds |
| P8.4 | Warp energy requirements | VN | Pfenning and Ford 1997; Van Den Broeck 1999 | Numbers change, sign does not |
| P9.1 | Horizon divergence | VN (1+1) | Finazzi et al. 2009 | Stability requirement only |
| P10.1 | Casimir, squeezing | PF (experiment) | Lamoreaux 1997; Vahlbruch et al. 2016 | None |
| P11 | G/c⁴ | PG | CODATA | None |
| P12.1 | Wormhole flare-out | PG | Hochberg and Visser 1997 | None |

The chain has exactly one premise that an FTL spacetime could satisfy while leaving every proven result intact: P6.3 on a self-consistent semiclassical background. P2.3, P6.1 and P8.3 for interacting fields are also unproven, but each of them, if false, makes FTL harder to realize or leaves a proven link in place.

## 3. Independent check of every metric shipped with Warp Factory

All five warp generators behave as published at their shipped example parameters: Alcubierre, Van Den Broeck, Modified Time and Lentz violate the null, weak, strong and dominant conditions, and the Fuchs et al. (2024) shell satisfies all four at its example velocity of 0.02 c. With its mass held fixed, the shell fails the null condition from about 0.03 c upward; Fuchs et al. state that such a threshold exists and leave its value to future work; the value is the one number below that is not in the published record.

**Method.** Warp Factory (Helmerich et al. 2024) is MATLAB. Its repository (NerdsWithAttitudes/WarpFactory, main branch, cloned Oct 20, 2026, MIT licence) ships the generators Minkowski, Schwarzschild, Alcubierre, Van Den Broeck, Modified Time, Lentz and Warp Shell. I re-implemented its pipeline in Python with the warp generators ported line by line, including the Lentz piecewise template and the Warp Shell TOV construction with its moving-average smoothing. Each metric is evaluated in its Galilean comoving frame on the z = 0 plane of a five-plane slab (one time slice, as the comoving generators require), with 4th-order central differences, T\_μν = (c⁴/8πG) G\_μν in SI with x⁰ = ct, an Eulerian orthonormal tetrad, 400 null directions and 10 rapidities up to 2.0 for the sampled conditions (Warp Factory defaults: 100 and 10), plus an exact Hawking-Ellis eigenvalue test that Warp Factory does not perform. Minima are taken over the region where |T| exceeds 10⁻⁸ of its maximum, excluding a four-cell margin. Validation is in Appendix B.

| Generator, shipped example parameters | ρ\_E min \[J/m³\] | ∫ρ\_E dV \[J\] | NEC min, sampled \[J/m³\] | Verdict | Published claim |
| --- | --- | --- | --- | --- | --- |
| Alcubierre comoving, v = 0.5 c, R = 5 m, σ = 0.5 m⁻¹, h = 0.1 m | −1.93×10⁴⁰ (max 3×10²⁵, roundoff) | −1.134×10⁴³ (closed form −1.134×10⁴³) | −3.9×10⁴¹ | all four violated | violated (Alcubierre 1994); match |
| Alcubierre, v = 0.9 c | −6.26×10⁴⁰ | −3.67×10⁴³ | −7.7×10⁴¹ | all four violated | match |
| Alcubierre, v = 2.0 c | −3.09×10⁴¹ | −1.81×10⁴⁴ | −2.2×10⁴² | all four violated | match |
| Van Den Broeck comoving, v = 0.1 c, R1 = 2 m, σ1 = 1 m⁻¹, R2 = 5 m, σ2 = 1 m⁻¹, A = 0.5, h = 0.05 m | −7.08×10⁴¹ (max +1.79×10⁴² in the B pocket) | −6.19×10⁴³ | −9.6×10⁴¹ | all four violated | violated (Van Den Broeck 1999); match |
| Modified Time comoving, v = 0.5 c, R = 5 m, σ = 0.5 m⁻¹, A = 2 (lapse 1/2 inside) | −3.68×10⁴⁰ | −1.91×10⁴³ | −4.7×10⁴¹ | all four violated | violated (Warp Factory paper); match |
| Modified Time, A = 0.5 (lapse 2 inside) | −9.02×10³⁹ | −6.11×10⁴² | −6.4×10⁴¹ | all four violated | match |
| Lentz comoving, v = 0.1 c, template scale 30/7 m, raw piecewise-constant shift as shipped | −3.69×10⁴² (max +9.8×10⁴², grid-regularized surface terms) | −1.36×10⁴¹ per metre of z | −8.3×10⁴³ | all four violated | Lentz 2021 claimed positive energy; Santiago, Schuster and Visser 2022 and the Warp Factory paper: violated. Match with the latter |
| Lentz, shift smoothed with a 0.3 m Gaussian | −8.23×10⁴⁰ (max +2.9×10⁴¹) | −2.2×10³⁸ per metre of z | −8.1×10⁴¹ | all four violated | same |
| Warp Shell comoving (Fuchs et al. 2024), R1 = 10 m, R2 = 20 m, M = 4.49×10²⁷ kg (2GM/R2c² = 1/3), smoothFactor 4000, h = 0.2 m, static (v = 0) | −1.7×10³⁴ against max +1.36×10⁴⁰ (ratio 10⁻⁶) | +4.11×10⁴⁴ over coordinate volume (Mc² = 4.03×10⁴⁴) | −2.5×10³⁴ (ratio −2×10⁻⁶) | all four hold within tolerance, type I everywhere | consistent |
| Warp Shell, v = 0.02 c (shipped example) | −1.7×10³⁴ | +4.11×10⁴⁴ | −2.5×10³⁴ | all four hold within tolerance | satisfied (Fuchs et al. 2024); match |
| Warp Shell, v = 0.2 c, mass unchanged | −1.7×10³⁴ (ρ\_E stays positive) | +3.95×10⁴⁴ | −8.7×10⁴⁰ (ratio −2.15 of max abs T = 4.0×10⁴⁰) | NEC, WEC, SEC, DEC violated over 31% of the matter region | threshold anticipated by Fuchs et al., value left to future work; see below |
| Warp Shell, v = 0.5 c | −7.2×10³⁸ | +3.16×10⁴⁴ | −2.7×10⁴¹ (ratio −2.9) | violated over 37% | see below |
| Warp Shell, v = 0.9 c | −2.6×10⁴⁰ | +1.49×10⁴⁴ | −4.7×10⁴¹ | violated | see below |

For the Alcubierre rows the Eulerian density is non-positive everywhere (the maxima are roundoff at 10⁻¹⁵ of the scale), and the integrated energy matches the closed form −(v²/12)(c⁴/G)∫f'² r² dr to 0.05%. The totals are large: −1.1×10⁴³ J at v = 0.5 c for a 5 m bubble is about −20 Earth masses of negative energy. For Van Den Broeck the B pocket adds a positive region, and the net is still negative. The Lentz figures are grid-dependent because the shipped shift is piecewise constant, so its stress tensor is a sum of surface distributions; the smoothed row is the meaningful one and the sign is the same in both.

### The Warp Shell velocity threshold

Same shell, same mass, same smoothing; only the shift amplitude v changes unless noted. NEC min is the sampled minimum over the matter region divided by max|T|; the last column is the fraction of matter-region points where the null condition fails by more than 10⁻³ max|T|.

| v \[c\] | 2GM/R2c² | Shift buffer R\_buff \[m\] | NEC min / max abs T | Fraction of matter region violating |
| --- | --- | --- | --- | --- |
| 0.02 | 0.33 | 0 | −2×10⁻⁶ | 0.000 |
| 0.05 | 0.33 | 0 | −0.64 | 0.062 |
| 0.10 | 0.33 | 0 | −1.61 | 0.190 |
| 0.20 | 0.33 | 0 | −2.15 | 0.317 |
| 0.20 | 0.50 | 0 | −2.12 | 0.276 |
| 0.20 | 0.70 | 0 | −2.23 | 0.240 |
| 0.20 | 0.33 | 3 | −2.32 | 0.218 |
| 0.50 | 0.70 | 0 | −2.26 | 0.332 |

At v = 0.2 c the worst point sits in the equatorial plane at r = 12.4 m. Its Eulerian-frame stress tensor has T₀₀ = 1.13×10⁴⁰, T₀₁ = −4.03×10⁴⁰, T₁₁ = −1.76×10⁴⁰, T₃₃ = −1.43×10⁴⁰ J/m³; the eigenvalues of T^a\_b are (−1.44 ± 4.02 i)×10⁴⁰, −0.04×10⁴⁰ and −1.43×10⁴⁰, so the tensor is Hawking-Ellis type IV and has no rest frame. Halving the grid spacing to 0.1 m changes none of these numbers to three digits. Violating points fill the wall from r = 10.1 m to 19.8 m.

My reading, not a theorem: in the wall the shift gradient contributes a momentum density linear in v (through the momentum constraint 8π j\_i = D\_j(K^j\_i − δ^j\_i K) with K\_ij ∝ v) and an energy density quadratic in v, while the shell's energy density does not depend on v. Along the null vectors k = n ± e the condition compares ρ\_E + p(e,e) with 2|j·e|, so it fails once the shift-induced momentum density reaches about half the shell density. For this shell that happens near 0.03 c. Raising the compactness to 0.5 and 0.7 raised the shell density and the shift-induced momentum together, and the threshold did not move above 0.2 c in those two trials; confining the shift change to the inner 4 m of the wall made it worse. Whether another radius ratio or shift profile moves the threshold is a Phase 1 calculation. Fuchs et al. (2024) report that the shift was kept below the level at which the added momentum flux causes violations, that 0.02 c works, and that the upper limit is future work; they also give the same momentum-against-density comparison as a rule of thumb. The new content here is the measured value of the threshold for the shipped shell and the observation that compactness alone did not move it.

### Repository observations (main branch, Oct 20, 2026)

- `metricGet_VanDenBroeckComoving` builds g\_tt from the shape function f where the Galilean comoving form of Van Den Broeck's line element needs (1 − f): g\_tt = −1 + B²v²(1 − f)². The difference is O(v²); at v = 0.1 c it changes the energy-condition minima by under 1% and no verdict. Both forms are in the table's source data.
- `metricGet_ModifiedTimeComoving` has the same g\_tt mismatch and two indexing slips: x is built from the time-axis grid scale and centre (`gridScaling(1)`, `worldCenter(1)`), and g₂₁ is copied from element (t, i, k, k) instead of (t, i, j, k). With A = 2 and v = 0.5 c the as-coded metric is degenerate inside the bubble (lapse² = 1/A² − v² = 0). The table uses the corrected comoving form; the as-coded A = 0.5 case gives the same verdict (NEC min −5.8×10⁴¹).
- The non-comoving generators are self-consistent. None of this affects any conclusion of the Warp Factory paper that I checked.

## 4. The exact loophole

Every proven link in Section 2 survives an FTL spacetime only if one premise fails: the averaged null energy condition along complete achronal null geodesics of a self-consistent semiclassical solution. Stated as a premise:

**P6.3 (self-consistent achronal ANEC).** Let (M, g) be a globally hyperbolic, asymptotically flat spacetime solving G\_μν = 8π ⟨T\_μν⟩\_ω, where ω is a state of the full interacting quantum field theory of the matter present and ⟨T\_μν⟩\_ω is renormalized by the Hadamard prescription. Then for every complete achronal null geodesic γ with affine parameter λ and tangent k,

```latex
\int_{\gamma} \langle T_{\mu\nu} \rangle_{\omega}\, k^{\mu} k^{\nu}\, d\lambda \;\geq\; 0 .
```

**Why this premise and no other.** Olum's theorem (L4) is unconditional in classical GR: a superluminal path needs pointwise NEC violation somewhere along it. Quantum fields can violate the NEC pointwise, so L4 alone does not close the door; it tells you where the negative energy must sit. Gao and Wald (L6) close the door with the integrated condition, and the integrated condition is proven only in flat space. The Natário class is closed by L7 whatever the state, so a candidate needs a non-unit lapse or non-flat slices, and Section 3 shows the published constructions of that kind satisfy the classical conditions only below a few percent of c with planetary masses. A candidate therefore has to be a self-consistent semiclassical solution whose negative energy survives the null average along the fastest achronal geodesic. That is a violation of P6.3 and of nothing else that is proven.

**Status of P6.3 in one table.**

| Setting | Result | Reference | Tag |
| --- | --- | --- | --- |
| Minkowski space, any unitary Lorentz-invariant QFT | Proven | Faulkner, Leigh, Parrikar and Wang 2016; Hartman, Kundu and Tajdini 2017 | PF |
| Free fields, weak curvature, classical NEC-satisfying background | Proven | Kontou and Olum 2015 | PG (restricted) |
| Semiclassical gravity assuming the generalized second law | Derived | Wall 2010 | CJ (GSL unproven) |
| Conformally flat background, not self-consistent, conformally coupled field | Violated along achronal geodesics | Urban and Olum 2010 | counterexample, outside the premise |
| Non-achronal geodesics (Casimir plates, compactified directions) | Violated | Graham and Olum 2007 and references therein | outside the premise |
| General self-consistent interacting case | Open | Graham and Olum 2007 conjecture | CJ |

**If P6.3 is proven.** L6 and topological censorship apply to every globally hyperbolic semiclassical spacetime: no null geodesic reaches infinity ahead of its flat-space counterpart and no traversable wormhole shortcut exists. FTL by spacetime shortcut is then excluded wherever semiclassical gravity is valid. What remains is full quantum gravity, where neither g nor ⟨T⟩ is a valid description, and non-globally-hyperbolic spacetimes, which carry Cauchy horizons or closed timelike curves and return the question to P2.3.

**If P6.3 is false.** A counterexample is a background on which a superluminal arrival is not forbidden by any theorem in Section 2. It is still only a candidate. To become a drive it must also satisfy the quantum inequalities pointwise (P8.3), which confine negative energy to thin regions and set the Pfenning and Ford and Van Den Broeck scale of the total unless the geometry is specifically built to evade it; avoid the horizon divergence of L9 if v exceeds c; survive P2.3 once two devices exist; and be a state that can be prepared, since a mathematically allowed state is not a laboratory procedure. None of these come with the counterexample.

**The other unproven premises, and why they are not loopholes.** P2.3 (chronology protection): if false, FTL devices are time machines, which does not make them easier to build. P6.1 (global hyperbolicity): if false, the spacetime has a Cauchy horizon or closed timelike curves and inherits P2.3. P8.3 for interacting fields: if the quantum inequalities failed for interacting matter, local negative energy would be less constrained, but L4 and L6 would stand, and the flat-space ANEC and QNEC results already hold for interacting theories, which points the other way.

**Gate A criteria, restated against P6.3.** A Track B candidate passes only if it exhibits a self-consistent semiclassical solution with achronal ANEC violation along a complete null geodesic that reaches infinity, satisfies every proven quantum inequality, carries no horizon, has total negative energy below about 1 J, and comes with a written matter Lagrangian and state. The first clause is the new one; the other four were in the Phase 0 plan.

## 5. Referee notes and Phase 1 targets

The no-go is complete at the classical level (L4, L7, L12) and at the flat-space quantum level (L6 with P6.3 in Minkowski space, L8). It is open at exactly one point, P6.3 on curved self-consistent backgrounds, and that point is where Phase 1 money goes.

**Settled by this document.** The twelve-link chain with every premise tagged (Section 2). Independent confirmation that all four Alcubierre-family generators in Warp Factory violate every pointwise energy condition at their shipped parameters, with totals that match closed forms to 0.05% (Section 3). Confirmation that the Fuchs et al. shell satisfies all four at 0.02 c. Two identities for the Natário class, verified to machine precision pointwise and to 10⁻⁴ integrated (Appendix A); both are likely implicit in Santiago, Schuster and Visser 2022 and are not claimed as new.

**Open after this document.** P6.3 in general (CJ). P2.3 chronology protection (CJ). P6.1 global hyperbolicity (CJ). P8.3 for interacting four-dimensional fields (CJ). The dependence of the Warp Shell velocity threshold on radius ratio and shift profile, which Fuchs et al. list as future work.

**Points a referee will raise, with my answers.**

1. The Lentz template is piecewise constant, so its stress tensor is distributional and the raw numbers depend on the grid. Agreed; the table reports the smoothed version beside it and the sign of the net energy is the same in both, which is what L7 predicts for any non-gradient shift.
2. The Warp Shell profiles are smoothed after the TOV step, so the static shell is not an exact TOV solution. Agreed; the residual energy-condition margins at v ≤ 0.02 c are 10⁻⁶ of the scale, consistent with zero at this resolution, and the type I classification holds at every point.
3. Two comoving generators in the repository carry a g\_tt mismatch. Reported in Section 3 with both forms evaluated; no verdict changes.
4. The Warp Shell threshold near 0.03 c is a numerical observation plus an interpretation. The observation converged under grid refinement; the interpretation (linear-in-v momentum against v-independent shell density) coincides with the rule of thumb Fuchs et al. state, and turning it into a bound is the first Track B calculation below.
5. The tag PG for P6.2 (null generic condition) is generous. It fails on exactly symmetric geometries; an FTL candidate built on a pp-wave or similar would need a separate check of L6.
6. Experimental tags (P1.1, P10.1) are PF because the analyses assume a flat background; this is the correct tag for what those experiments test, and nothing in the chain needs more.

**Phase 1 targets, Jan 20, 2027 to Oct 20, 2027.**

- [ ] Track A1: prove P6.3 for free fields on self-consistent backgrounds beyond perturbation theory, extending Kontou and Olum 2015 to backgrounds that themselves carry quantum stress-energy. A proof here ends the semiclassical question for free matter.
- [ ] Track A2: attempt a self-consistent counterexample built on the Urban and Olum 2010 mechanism, with backreaction included and the geodesic checked for achronality.
- [ ] Track B1: derive the mass-velocity relation for the Fuchs et al. shell analytically from the momentum constraint, then test it against the Section 3 scan at three radius ratios. This settles whether the 0.03 c threshold is a property of the construction or of the parameters.
- [ ] Track B2: apply Fewster's curved-space quantum inequality (P8.3, free fields) to the Modified Time and Warp Shell wall regions where the null condition fails above threshold, to bound how much of the violation a quantum state could supply.
- [ ] Track B3: candidate search restricted to non-unit-lapse, non-flat-slice geometries, since L7 closes the Natário class for every state; each candidate runs through the Section 3 pipeline and the Gate A criteria of Section 4.
- [ ] Done: Fuchs et al. 2024 report 0.02 c as a working value, state that an upper limit on the shift exists, and leave it to future work. No discrepancy with Section 3.

**Gate A, Oct 20, 2027.** Proceed only with a candidate meeting all five criteria at the end of Section 4. My expectation from this document is that the gate fails and the program publishes the no-go with the P6.3 status table as its central claim.

## Appendix A. Eulerian energy of the Natário class

For unit lapse and flat slices the Eulerian energy density integrates to a non-positive total, and the direction-averaged null energy integrates to exactly twice that total. Both statements are derived here and both were checked on the Section 3 data.

**Setup.** Metric ds² = −dt² + δ\_ij (dx^i − β^i dt)(dx^j − β^j dt), G = c = 1. The slices t = const are flat, so the extrinsic curvature and its trace are

```latex
K_{ij} = \tfrac{1}{2}\left(\partial_i \beta_j + \partial_j \beta_i\right), \qquad K = \partial_i \beta^i .
```

The Hamiltonian constraint with ³R = 0 gives the Eulerian energy density

```latex
16\pi\, \rho_E = K^2 - K_{ij} K^{ij} .
```

**Step 1, the total.** Expanding K\_ij K^ij = ½(∂\_iβ\_j ∂\_iβ\_j + ∂\_iβ\_j ∂\_jβ\_i) and integrating the second term by parts twice (the derivatives of β must decay at infinity; β itself may tend to a constant, as it does in the comoving frame) gives ∫ ∂\_iβ\_j ∂\_jβ\_i = ∫ K², hence

```latex
16\pi \int \rho_E\, d^3x = \tfrac{1}{2} \int \left( K^2 - \partial_i \beta_j\, \partial_i \beta_j \right) d^3x .
```

**Step 2, the sign.** In Fourier space K̂ = i k·β̂ and |k·β̂|² ≤ |k|² |β̂|² by Cauchy-Schwarz, so ∫K² ≤ ∫ ∂\_iβ\_j ∂\_iβ\_j with equality only when β̂ is parallel to k for almost every k, that is, when β is a gradient. Therefore ∫ρ\_E d³x ≤ 0, with equality only for β = ∇φ.

**Step 3, the gradient case.** For β = ∇φ the density becomes

```latex
16\pi\, \rho_E = (\nabla^2 \phi)^2 - (\partial_i \partial_j \phi)(\partial_i \partial_j \phi) = 2 \sum_{i<j} \lambda_i \lambda_j ,
```

with λ\_i the eigenvalues of the Hessian of φ. It integrates to zero and vanishes identically only if the Hessian has rank at most one everywhere, which forces φ = h(a·x): a plane-symmetric shift that is not localized. So every localized non-trivial shift in the class has ρ\_E < 0 on an open set, and the weak and dominant conditions fail for Eulerian observers.

**Step 4, the null average.** Let S be the trace of the Eulerian spatial stress. For a stationary metric with unit lapse the trace of the ADM evolution equation reads β^i ∂\_i K = −K\_ij K^ij − 4π(ρ\_E + S). Combining with the constraint,

```latex
\rho_E + \tfrac{1}{3} S = \frac{1}{24\pi}\left( K^2 - 3 K_{ij}K^{ij} - 2\, \beta^i \partial_i K \right).
```

Integrating with ∫ β^i ∂\_i K = −∫ K² gives

```latex
\int \left(\rho_E + \tfrac{1}{3} S\right) d^3x = \frac{1}{8\pi} \int \left(K^2 - K_{ij}K^{ij}\right) d^3x = 2 \int \rho_E\, d^3x \;\leq\; 0 .
```

The left side is the average over spatial unit vectors e of T\_μν k^μ k^ν for the null vectors k = n + e, because the flux term 2 j·e averages to zero. For any non-gradient shift the averaged null energy is therefore negative on an open set, and the null energy condition fails there. Santiago, Schuster and Visser (2022) prove the NEC statement for the whole class including the gradient case; the derivation above is the short route for the non-gradient case and is what Section 3's generators all exercise.

**Numerical check.** On the Alcubierre data of Section 3 (v = 0.5 c, R = 5 m), the pipeline's ρ\_E matches the constraint formula to 3×10⁻¹⁵ and its S matches the Step 4 identity to 1.8×10⁻⁵ pointwise. The integrated ratio ∫(ρ\_E + S/3) / ∫ρ\_E is 2.0001 for a wall with σ = 1.5 m⁻¹ on a ±10 m domain, and 2.084 for σ = 0.5 m⁻¹ on the same domain, where the truncated tails account for the difference.

**What this appendix does not cover.** Non-unit lapse (Modified Time, Warp Shell), non-flat slices (Van Den Broeck), and time-dependent shifts in Step 4; Steps 1 to 3 hold at each instant for a time-dependent β on flat slices. Those cases are handled numerically in Section 3.

## Appendix B. Numerical method, validation, reproducibility

The pipeline reproduces an exact symbolic Einstein tensor to 1.5%, the Alcubierre closed-form energy density to 7×10⁻⁶, and a vacuum solution to 4×10⁻⁴ of its curvature scale; those three numbers bound the trust to place in Section 3.

**Pipeline.** Each generator returns the covariant metric on a slab of Nx × Ny × 5 points, spacing h in all three directions, in the Galilean comoving frame (one time slice; time derivatives vanish). First and second metric derivatives at the central plane use 4th-order central stencils (coefficients ±1, ∓8 over 12h for first derivatives; −1, 16, −30, 16, −1 over 12h² for second; mixed derivatives by composing first-derivative stencils). The Ricci tensor uses the coordinate formula R\_μν = ½ g^αβ(∂\_α∂\_ν g\_μβ + ∂\_μ∂\_β g\_αν − ∂\_α∂\_β g\_μν − ∂\_μ∂\_ν g\_αβ) + g^αβ g^γδ(Γ\_αμγ Γ\_βνδ − Γ\_αμν Γ\_βγδ), with the first index of Γ lowered; the ΓΓ term was checked symbolically against the Christoffel definition on a generic four-function metric before use (a first version had the wrong second term and was caught by that check). Then G\_μν = R\_μν − ½R g\_μν and T\_μν = (c⁴/8πG) G\_μν with c = 299 792 458 m/s and G = 6.674 30×10⁻¹¹ m³ kg⁻¹ s⁻², the values Warp Factory uses. The Eulerian tetrad is e₀ = (1/α, −β^i/α) with spatial legs from the Cholesky factor of γ\_ij. Energy conditions are evaluated in that frame on 400 Fibonacci-distributed null directions and 10 rapidities from 0 to 2.0 (null: T(k,k); weak and strong: T(u,u) and T(u,u) + ½T; dominant: the vector −T^a\_b k^b checked future-causal). The eigenvalue test diagonalizes T^a\_b = η^ac T\_cb, identifies the single timelike eigenvector where one exists (Hawking-Ellis type I), and reads the conditions off ρ and the three principal pressures: NEC ρ + p\_i ≥ 0, WEC adds ρ ≥ 0, DEC ρ ≥ |p\_i|, SEC adds ρ + Σp\_i ≥ 0. Points with complex eigenvalues or no timelike eigenvector are reported as non-type-I; the sampled test covers them.

**Validation.**

| Test | Result |
| --- | --- |
| Minkowski slab | every component of G exactly 0 |
| Schwarzschild, isotropic coordinates, r\_s = 0.2 m, h = 0.1 m, r > 1 m | max abs G = 7.3×10⁻⁵ against a curvature scale r\_s/r³ = 0.2 at r = 1 m (4×10⁻⁴ relative) |
| Alcubierre, v = 0.5 c, R = 5 m, σ = 0.5 m⁻¹, h = 0.1 m: full G\_μν against an exact sympy Einstein tensor computed by the Christoffel route | max abs difference 5.5×10⁻⁴ on a scale of 3.7×10⁻² (1.5%, finite-difference truncation in the wall) |
| Same: Eulerian ρ\_E against −v² y² f'² / (32π r²) × c⁴/G | 6.6×10⁻⁶ relative |
| Same: ρ\_E against the Hamiltonian constraint built from the same shift | 3×10⁻¹⁵ relative |
| Same: total Eulerian energy against −(v²/12)(c⁴/G) ∫ f'² r² dr | −1.1337×10⁴³ J against −1.1343×10⁴³ J (0.05%) |
| Appendix A identities | S\_E to 1.8×10⁻⁵ pointwise; integrated ratio 2.0001 |
| Warp Shell, v = 0.2 c, h = 0.2 m against 0.1 m | NEC minimum, location and frame tensor unchanged to three digits |

**Conventions that matter when comparing with Warp Factory output.** Warp Factory's `getEnergyTensor` returns contravariant T^μν in the coordinate basis; the Section 3 numbers are frame-tensor components T\_(a)(b) in the Eulerian orthonormal frame, which is what `doFrameTransfer(…, "Eulerian")` produces before `getEnergyConditions` samples them. Signs of the Eulerian momentum density follow T\_(0)(i) directly; Warp Factory's plotting scripts flip that sign for display. Edge cells are excluded here (four-cell margin) rather than filled by copying, which Warp Factory's `takeFiniteDifference1` does.

**Reproducibility.** Six scripts, delivered alongside this document: `wfcheck.py` (generators, derivatives, curvature, tetrad, energy conditions), `validate.py` (the symbolic comparison and vacuum tests), `run_checks.py` (every row of the main table; writes `results.json`), `shell_diag.py` and `shell_scan.py` (the Warp Shell threshold study), `appendixA_check.py` (the identities), plus `formula_test2.py` (the symbolic check of the Ricci formula). Runtime for `run_checks.py` is about three minutes on one CPU core with 3 GB of memory; the Warp Shell cases dominate. Parameters are those of the Warp Factory example live scripts: Alcubierre and Modified Time from `M2_Default_Metrics.mlx`, Van Den Broeck and Lentz likewise, Warp Shell from `W1_Warp_Shell.mlx` (R1 = 10, R2 = 20, m = R2 c²/(2G) × 1/3, vWarp = 0.02, sigma = 0, smoothFactor = 4000, spaceScale = 5). Where the example used a coarser grid (one grid point per metre for Alcubierre and Van Den Broeck), the spacing here is 0.1 m and 0.05 m; the shipped spacing of 0.2 m is kept for the Warp Shell.

## Sources

Author, year, title and journal are given; volume and page numbers are omitted so that nothing uncertain is printed as fact. Code: [Warp Factory repository](https://github.com/NerdsWithAttitudes/WarpFactory), main branch, cloned Oct 20, 2026.

- Abdo A.A. et al. (Fermi LAT and GBM) 2009, A limit on the variation of the speed of light arising from quantum gravity effects, Nature.
- Alcubierre M. 1994, The warp drive: hyper-fast travel within general relativity, Class. Quantum Grav.
- Antonello M. et al. (ICARUS) 2012, Measurement of the neutrino velocity with the ICARUS detector at the CNGS beam, Phys. Lett. B.
- Balakrishnan S., Faulkner T., Khandker Z.U., Wang H. 2019, A general proof of the quantum null energy condition, JHEP.
- Benford G.A., Book D.L., Newcomb W.A. 1970, The tachyonic antitelephone, Phys. Rev. D.
- Bobrick A., Martire G. 2021, Introducing physical warp drives, Class. Quantum Grav.
- Bousso R., Fisher Z., Koeller J., Leichenauer S., Wall A.C. 2016, Proof of the quantum null energy condition, Phys. Rev. D.
- Everett A.E. 1996, Warp drive and causality, Phys. Rev. D.
- Everett A.E., Roman T.A. 1997, Superluminal subway: the Krasnikov tube, Phys. Rev. D.
- Faulkner T., Leigh R.G., Parrikar O., Wang H. 2016, Modular Hamiltonians for deformed half-spaces and the averaged null energy condition, JHEP.
- Fell S.D.B., Heisenberg L. 2021, Positive energy warp drive from hidden geometric structures, Class. Quantum Grav.
- Fewster C.J. 2000, A general worldline quantum inequality, Class. Quantum Grav.
- Fewster C.J., Eveson S.P. 1998, Bounds on negative energy densities in flat spacetime, Phys. Rev. D.
- Fewster C.J., Hollands S. 2005, Quantum energy inequalities in two-dimensional conformal field theory, Rev. Math. Phys.
- Finazzi S., Liberati S., Barceló C. 2009, Semiclassical instability of dynamical warp drives, Phys. Rev. D.
- Ford L.H., Roman T.A. 1995, Averaged energy conditions and quantum inequalities, Phys. Rev. D.
- Ford L.H., Roman T.A. 1997, Restrictions on negative energy density in flat spacetime, Phys. Rev. D.
- Friedman J.L., Schleich K., Witt D.M. 1993, Topological censorship, Phys. Rev. Lett.
- Fuchs J., Helmerich C., Bobrick A., Sellers L., Melcher B., Martire G. 2024, Constant velocity physical warp drive solution, Class. Quantum Grav.
- Gao S., Wald R.M. 2000, Theorems on gravitational time delay and related issues, Class. Quantum Grav.
- Geroch R. 1967, Topology in general relativity, J. Math. Phys.
- Graham N., Olum K.D. 2007, Achronal averaged null energy condition, Phys. Rev. D.
- Hartman T., Kundu S., Tajdini A. 2017, Averaged null energy condition from causality, JHEP.
- Hawking S.W. 1992, Chronology protection conjecture, Phys. Rev. D.
- Helmerich C., Fuchs J., Bobrick A., Sellers L., Melcher B., Martire G. 2024, Analyzing warp drive spacetimes with Warp Factory, Class. Quantum Grav.
- Hiscock W.A. 1997, Quantum effects in the Alcubierre warp-drive spacetime, Class. Quantum Grav.
- Hochberg D., Visser M. 1997, Geometric structure of the generic static traversable wormhole throat, Phys. Rev. D.
- Hochberg D., Visser M. 1998, Dynamic wormholes, antitrapped surfaces, and energy conditions, Phys. Rev. D.
- Kay B.S., Radzikowski M.J., Wald R.M. 1997, Quantum field theory on spacetimes with a compactly generated Cauchy horizon, Commun. Math. Phys.
- Klinkhammer G. 1991, Averaged energy conditions for free scalar fields in flat spacetime, Phys. Rev. D.
- Kontou E.-A., Olum K.D. 2015, Proof of the averaged null energy condition in a classical curved spacetime using a null-projected quantum inequality, Phys. Rev. D.
- Krasnikov S.V. 1998, Hyperfast travel in general relativity, Phys. Rev. D.
- Lamoreaux S.K. 1997, Demonstration of the Casimir force in the 0.6 to 6 μm range, Phys. Rev. Lett.
- Lentz E.W. 2021, Breaking the warp barrier: hyper-fast solitons in Einstein-Maxwell-plasma theory, Class. Quantum Grav.
- Lobo F.S.N., Visser M. 2004, Fundamental limitations on warp drive spacetimes, Class. Quantum Grav.
- McMonigal B., Lewis G.F., O'Byrne P. 2012, Alcubierre warp drive: on the matter of matter, Phys. Rev. D.
- Mohideen U., Roy A. 1998, Precision measurement of the Casimir force from 0.1 to 0.9 μm, Phys. Rev. Lett.
- Morris M.S., Thorne K.S. 1988, Wormholes in spacetime and their use for interstellar travel: a tool for teaching general relativity, Am. J. Phys.
- Morris M.S., Thorne K.S., Yurtsever U. 1988, Wormholes, time machines, and the weak energy condition, Phys. Rev. Lett.
- Nagel M. et al. 2015, Direct terrestrial test of Lorentz symmetry in electrodynamics to 10⁻¹⁸, Nature Communications.
- Natário J. 2002, Warp drive with zero expansion, Class. Quantum Grav.
- Olum K.D. 1998, Superluminal travel requires negative energies, Phys. Rev. Lett.
- Pfenning M.J., Ford L.H. 1997, The unphysical nature of warp drive, Class. Quantum Grav.
- Santiago J., Schuster S., Visser M. 2022, Generic warp drives violate the null energy condition, Phys. Rev. D.
- Tipler F.J. 1977, Singularities and causality violation, Ann. Phys.
- Tolman R.C. 1917, The Theory of the Relativity of Motion, University of California Press.
- Urban D., Olum K.D. 2010, Averaged null energy condition violation in a conformally flat spacetime, Phys. Rev. D.
- Vahlbruch H., Mehmet M., Danzmann K., Schnabel R. 2016, Detection of 15 dB squeezed states of light and their application for the absolute calibration of photoelectric quantum efficiency, Phys. Rev. Lett.
- Van Den Broeck C. 1999, A warp drive with more reasonable total energy requirements, Class. Quantum Grav.
- Visser M., Bassett B., Liberati S. 2000, Superluminal censorship, Nucl. Phys. B Proc. Suppl.
- Wall A.C. 2010, Proving the achronal averaged null energy condition from the generalized second law, Phys. Rev. D.
