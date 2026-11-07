# Ricci Soliton Rigidity

Research on complete steady gradient Ricci solitons with the standard
\(O(p)\times O(q)\) symmetry, especially the four-dimensional
\(p=q=2\) Lavoyer--Peachey family.

## Established in the working proof

* **Boundary-orbit criterion:** noncollapse is equivalent to positive
  asymptotic orbit-radius slope on the fixed submanifold. For \(p=2\),
  \(\alpha>0\) is equivalent to total Gaussian curvature of that fixed
  surface being strictly less than \(2\pi\). The proof compares intrinsic
  and ambient geometry without assuming global convexity of the fixed set.
* **Exact symmetry in all factor dimensions:** distinct pure sectional
  curvatures at the tip determine the full isometry group through the
  curvature operator, even when the Ricci eigenvalues coincide. Thus
  \(\operatorname{Isom}(\mathcal S^{p,q}_\theta)=O(p)\times O(q)\) for
  all \(p,q\ge2\), \(0<\theta<1\).
* **Four-dimensional collapse corollary:** assuming Theorem 1.1 of Fang's
  27 October 2026 preprint, \(\mathcal S^{2,2}_\theta\) is not
  \(\kappa\)-noncollapsed at all scales for any \(\kappa>0\). The proof
  is an incompatible-isometry-group argument; the analytic symmetry
  theorem is an explicit external dependency.
* **Correct reduction:** the ambient action has cohomogeneity two, not one.
  The exact two-variable soliton equations and weighted identities are
  derived in the proof.
* **Flat-base rigidity:** under nonnegative sectional curvature, a flat
  orbit-space metric forces a product. In the positively Ricci-curved
  four-dimensional steady case the factors are scaled cigars. These
  products have arbitrary tip-curvature ratio but collapsed boundary
  orbits, and demonstrate why strict curvature cannot be weakened.
* **Convergent local-soliton obstruction:** explicit collapsed axis data
  extend to genuine analytic steady solitons on neighborhoods of a complete
  fixed surface, smoothly through the vertex, with positive curvature
  operator and any prescribed tip ratio. The surface has total curvature
  \(2\pi\). The proof combines an all-orders recurrence, Bianchi constraint
  propagation, and an analytic convergence argument. A uniform-width
  positive-curvature tube exists along the entire surface; its end
  converges exponentially to a flat cylinder times a cigar. The first
  transverse end mode grows exponentially without a sign obstruction in
  its first-circle curvature block.
* **Quantitative nonlinear continuation at infinity:** there are
  constants \(s_*,C_*>0\) such that the exact positive-curvature solution
  extends throughout \(s\ge s_*+C_*t\), \(t\ge0\).
  Thus finite-width thresholds satisfy \(s_T\le s_*+C_*T\).
  The proof combines fixed-factor analytic Cauchy restarts with exponential
  curvature lower bounds and a Schur-complement estimate. The complementary
  region and general second singular-orbit matching remain unresolved. A fixed \(s\) still permits only a finite normal interval.

* **Necessary global curvature budget:** any complete positively curved steady soliton
  extension of the constructed germ must have reflected-section curvature
  \(\int_S K_S\,dA\le4\arctan(3/\sqrt{4/\theta-1})\).
  Gauss--Bonnet on potential-level quarter-disks proves the bound.
  This section is not the fixed surface in the boundary-orbit criterion;
  the condition does not yet exclude complete extension.

* **Small-parameter obstruction to completion:** for sufficiently small
  \(\theta\), the constructed cigar-axis germ develops negative mixed-circle
  sectional curvature at \(t=2A\log(1/A)\), where
  \(A=2/\sqrt{4/\theta-1}\). A uniform logarithmic-time analytic estimate
  proves this is an actual curvature sign change. Consequently this germ
  has no complete positively curved steady soliton extension in that range.
  This does not identify the axis data of the Lavoyer--Peachey examples
  or prove their strict curvature gap.

See [PROOF.md](PROOF.md) for complete arguments, source statements, and
precise hypotheses. No numerical experiment is used.

## Still open

The main target remains proving the strict total-curvature bound, hence
\(\alpha=a(0)>0\), for Lavoyer--Peachey's examples. This is
**boundary-orbit noncollapse**, not ambient
\(\kappa\)-noncollapse. The closed-disk conclusion and the paper's enhanced
dimension reduction remain conditional on that hypothesis.

These notes make no priority claim and are not a submission manuscript.
