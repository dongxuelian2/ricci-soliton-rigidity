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
* **Actual-axis alternatives:** the axis radius is bounded exactly when
  its angular--normal curvature is integrable. An unbounded sublinear
  radius requires nonsummable angular--normal curvature but summable
  radial--normal curvature. Explicit analytic axis data satisfy all these
  tests and every tip ratio. Proposition 6.7 now excludes complete
  positive-curvature extensions of that entire diagnostic family by an
  ambient sectional-curvature sign obstruction, including the previously
  unexcluded parameter range. Thus those axis tests remain insufficient,
  but that family supplies no candidate for complete soliton existence.
* **Actual ambient limits and exclusion of analytic sublinear power laws:** every
  unbounded actual axis, including a sublinear one, has unrescaled
  pointed limits equal to a plane times a cigar, or flat four-space.
  These are actual products, with quotient possibilities excluded.
  If the limiting normal sectional curvature is positive, an analytic
  power-law expansion
  \(w=w_0s^\mu(1+\omega/s+O(s^{-2}))\),
  \(F_0'=\ell-b/s+g/s^2+O(s^{-3})\) with \(0<\mu<1\)
  is impossible. Two sectional-curvature tests first force exact
  equalities for \(b,g\). An exact adjoint integral evaluates the second
  transverse coefficient. The remaining meridional curvature is then
  strictly negative at order \(s^{-4}\), independently of all higher
  axis coefficients. A fourth-order development with a uniform error
  on \(t=\gamma\log s\) justifies the sign for the actual soliton.
  This replaces the former interval restriction by exclusion of the
  entire analytic power-law sublinear regime. Corollary 6.7a also permits
  finite fourth-order holomorphic asymptotics with a uniformly controlled
  remainder, without a convergent series at infinity. The LP expansion
  hypothesis remains unproved.
* **Global transverse curvature budget:** every normal geodesic from
  one fixed axis meets the other axis in finite time, with nonorthogonal
  terminal tangent. The first-circle Jacobi equation gives
  \(\int_0^T K a\,dt<1\) and
  \(\int_0^{w/2}K\,dt<2/w\), without curvature-scale compactness or
  a globally injective normal gauge. A separate, still unproved
  transverse lower bound would restrict power-law exponents to
  \(\mu\le1/2\); the axis curvature value alone does not supply it.
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
  the budget alone does not exclude complete extension.

* **Small-parameter obstruction to completion:** for sufficiently small
  \(\theta\), the constructed cigar-axis germ develops negative mixed-circle
  sectional curvature at \(t=2A\log(1/A)\), where
  \(A=2/\sqrt{4/\theta-1}\). A uniform logarithmic-time analytic estimate
  proves this is an actual curvature sign change. Consequently this germ
  has no complete positively curved steady soliton extension in that range.
  This does not identify the axis data of the Lavoyer--Peachey examples
  or prove their strict curvature gap.

* **Nonextension for every parameter:** the end-mode argument excludes
  complete positively curved steady soliton extensions of this same germ
  for every \(0<\theta<1\). Uniform shifted-mode semigroup bounds and a
  shrinking analytic-radius estimate control the nonlinear remainder
  where the mixed-circle curvature becomes negative. The actual
  Lavoyer--Peachey axis profiles are not identified with this germ,
  so their strict-gap problem remains open.

See [PROOF.md](PROOF.md) for complete arguments, source statements, and
precise hypotheses. No numerical experiment is used.

## Still open

The main target remains proving the strict total-curvature bound, hence
\(\alpha=a(0)>0\), for Lavoyer--Peachey's examples. This is
**boundary-orbit noncollapse**, not ambient
\(\kappa\)-noncollapse. The closed-disk conclusion and the paper's enhanced
dimension reduction remain conditional on that hypothesis.

Proposition 6.7 eliminates every analytic sublinear power-law axis
with positive limiting normal curvature, not only the explicit
diagnostic family. To turn that exclusion into a strict gap for LP,
one must control the actual axis data and their normal development
on distances growing with the axis parameter, without assuming an
analytic expansion at infinity.
The fixed-radius product limit alone supplies neither that expansion
nor such uniform control. The branch with zero limiting normal
curvature also remains separate.

The proof now gives a conditional Riccati criterion for the strict
slope gap, together with oscillatory end data showing that the required
relative derivative bounds do not follow from the existing axis tests
and bounded axis derivatives. Deriving the criterion's hypothesis from
the actual soliton equations remains open.

The exact fourth-order source, adjoint integral and full curvature
coefficient can be reconstructed with
`python scripts/check_sublinear_fourth_order.py` using SymPy. This checks
the algebra; it does not verify the analytic hypotheses for actual LP axes.

These notes make no priority claim and are not a submission manuscript.
