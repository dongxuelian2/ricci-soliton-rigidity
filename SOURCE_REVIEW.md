# Review of a product-limit alternative in the LP construction

## Scope

This note checks one step relevant to the original boundary-orbit question: whether the tip-curvature ratio selected in the Lavoyer–Peachey construction, together with smooth convergence and the stated symmetry, excludes a product of two cigars.

It does **not** prove that the actual LP sequence converges to a product, disprove the existence theorem, or supply a counterexample with strictly positive curvature operator. The strict total-curvature gap remains unresolved.

Source: [Lavoyer–Peachey, published version](https://doi.org/10.1007/s12220-026-02567-x), Section 4, Lemma 4.1 and the proof of Theorem 1.1; see also [arXiv:2509.25081v2](https://arxiv.org/html/2509.25081v2#S4).

## What the inspected argument supplies

The construction takes a locally smooth limit of tip-normalized expanding solitons and uses an intermediate-value argument to impose a prescribed ratio between the two pure sectional curvatures at the vertex. These operations preserve nonnegative curvature and the ratio. Strict positivity does not follow from smooth convergence alone.

No additional exclusion of the product alternative below was located in the inspected proof or its stated assumptions. This identifies an argument that still needs to be supplied or located, rather than asserting that the theorem is false.

## The remaining scalar degree of freedom

In dimension four, write the tangent representation at the vertex as
\[
T_oM=V_2\oplus W_2.
\]
Let \(x,y\) be the sectional curvatures of the two pure planes and \(z\) the common mixed sectional curvature. The full \(O(2)\times O(2)\) symmetry makes the curvature operator diagonal, with eigenvalues
\[
x,\quad y,\quad z,z,z,z.
\]
For example, reflections rule out the off-diagonal terms and rotations identify the mixed entries.

Under \(R(o)=1\) and \(x=\theta y\),
\[
2x+2y+8z=1,
\qquad
z=\frac{1-2(1+\theta)y}{8}.
\]
Thus prescribing \(0<\theta<1\) does not itself prevent \(z\) from tending to zero. The endpoint
\[
y=\frac1{2(1+\theta)},\qquad
x=\frac{\theta}{2(1+\theta)},\qquad z=0
\]
is compatible with these identities.

## An exact complete model realizing the endpoint

For \(k>0\), set
\[
w_k(s)=\sqrt{\frac2k}\tanh\!\left(\sqrt{\frac k2}\,s\right),
\quad
g_k=ds^2+w_k(s)^2d\phi^2,
\quad
F_k(s)=2\log\cosh\!\left(\sqrt{\frac k2}\,s\right).
\]
The usual smooth extension at \(s=0\) is a complete cigar. Direct differentiation gives
\[
K_k=-\frac{w_k''}{w_k}
=k\,\operatorname{sech}^2\!\left(\sqrt{\frac k2}\,s\right)
=F_k''
=F_k'\frac{w_k'}{w_k}.
\]
Hence \(\operatorname{Ric}(g_k)=\nabla^2F_k\).

Take
\[
a=\frac{\theta}{2(1+\theta)},\qquad
b=\frac1{2(1+\theta)}
\]
and the product \((g_a\oplus g_b,F_a+F_b)\). It has:

- completeness and the standard \(O(2)\times O(2)\) action;
- positive Ricci tensor and a unique potential critical point;
- nonnegative curvature operator, with four zero mixed eigenvalues;
- \(R(o)=2a+2b=1\) and the required pure curvature ratio \(a/b=\theta\).

Both fixed surfaces satisfy
\[
\int K\,dA
=2\pi\bigl(w_k'(0)-w_k'(\infty)\bigr)
=2\pi.
\]
This is why an argument using only those limiting properties cannot establish the desired strict bound. This product fails the theorem's strict curvature hypothesis.

## What must be checked on the actual sequence

At roots of the tip-ratio equation, ruling out this endpoint requires additional control of the mixed curvature. For instance, along the relevant convergent subsequence, a bound
\[
\limsup_m y_m<\frac1{2(1+\theta)}
\]
would imply \(\liminf_m z_m>0\). At a limit with \(x,y,z>0\), positivity at the vertex permits the usual strong maximum principle to propagate positivity. The maximum principle does not, by itself, eliminate a product when \(z=0\).

Such a mixed-curvature bound would address this product alternative; it would still not prove the strict axis total-curvature gap. No estimate of this kind for the actual selected sequence is established here.

The symmetry algebra, cigar identities, normalization and scope of this observation were checked separately. The two scalar steady equations and the total-curvature calculation also pass exact symbolic verification.

