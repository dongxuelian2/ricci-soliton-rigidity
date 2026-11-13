# Symmetry, collapse, and the boundary-orbit problem

## 1. Target and literature baseline

Write the steady equation as
\[
\operatorname{Ric}_g=\nabla^2 F.
\]
Thus the potential in the original repository convention is \(f=-F\).
All group actions below are equivariantly diffeomorphic to the standard
\(O(p)\times O(q)\) action on \(\mathbb R^p\oplus\mathbb R^q\).

Lavoyer--Peachey [LP, Theorem 1.1] construct positively curved examples
\(\mathcal S^{p,q}_\theta\), with \(p,q\ge2\), \(0<\theta<1\), and tip
sectional curvatures \(\lambda_p=\theta\lambda_q\).
Their Theorem 5.1 describes the ideal boundary as
\[
 [0,L]\times_a S^{p-1}\times_b S^{q-1},
 \qquad 0\le L\le\pi/2,
\]
with nonnegative concave Lipschitz warping functions, \(a\) decreasing,
\(b\) increasing, and \(a(L)=b(0)=0\).
Remark 5.2 leaves \(\alpha=a(0)>0\) open; if it holds, \(b\equiv0\) and
the boundary is a closed \(p\)-disk. The Euclidean splitting in Theorem 1.2
assumes this boundary conclusion. It cannot be used to prove its own
hypothesis. Remark 1.3 addresses ambient noncollapse only for \(q\ge3\).

Two different questions must therefore be kept separate:

* **Boundary-orbit noncollapse:** \(\alpha>0\), a distance-scale statement
  about the action on the ideal boundary.
* **All-scale curvature noncollapse:** existence of \(\kappa>0\) such that
  \[
  \sup_{B_g(x,r)}|\mathrm{Rm}|\le r^{-2}
  \quad\Longrightarrow\quad
  \operatorname{Vol}_g B_g(x,r)\ge\kappa r^{\dim M}.
  \tag{1}
  \]

The boundary question remains open in this document. Section 2 proves an
exact isometry-group statement without further literature input. Section 3
deduces failure of (1) in dimension four from Fang's recent symmetry theorem.
Sections 4--5 give the correct PDE and a rigidity obstruction to a flat-base
reduction. Section 6 proves an intrinsic criterion for the boundary question,
which becomes a strict total-curvature bound when \(p=2\).
These are reviewable deductions, not a claim to have resolved
\(\alpha>0\) or established research priority.

## 2. The full isometry group is determined at the vertex

### Proposition 2.1 (spectral determination of the isometry group)

Let \((M^{p+q},g,F)\) be a complete connected steady gradient Ricci soliton
with \(\operatorname{Ric}>0\), a critical point \(o\) of \(F\), and the standard
isometric \(O(p)\times O(q)\) action fixing \(o\). Write
\(T_oM=V_p\oplus V_q\). Let \(\lambda_p,\lambda_q\) denote the sectional
curvatures of planes in the two summands, and \(\mu\) the sectional curvature
of a plane spanned by one unit vector from each summand. Set
\[
 \rho_p=(p-1)\lambda_p+q\mu,
 \qquad \rho_q=(q-1)\lambda_q+p\mu.
 \tag{2}
\]
If \(\rho_p\ne\rho_q\), the given \(O(p)\times O(q)\) is the full isometry
group of \((M,g)\).

**Proof.** The contracted Bianchi identity and the steady equation give
\[
 \nabla R=-2\operatorname{Ric}(\nabla F,\cdot),
 \qquad R+|\nabla F|^2=R(o).
 \tag{3}
\]
Since \(\nabla^2F=\operatorname{Ric}>0\), \(F\) has at most one critical point:
along a minimizing geodesic between two distinct critical points,
\((F\circ\gamma)''>0\), contradicting zero derivatives at both endpoints.
Consequently \(R\) attains its global maximum at precisely \(o\), by (3).
Every isometry fixes \(o\), since it preserves scalar curvature.

The two summands are orthogonal and the Ricci endomorphism is scalar on each.
Indeed, the full orthogonal group on either summand acts transitively on its
unit sphere, and the element acting as \(-I\) on one summand and \(I\) on the
other eliminates cross terms. The same transitivity makes all mixed
sectional curvatures equal to \(\mu\). Tracing sectional curvatures gives
(2), so the distinct eigenspaces of the Ricci endomorphism are exactly
\(V_p,V_q\).

For any isometry \(\varphi\), its derivative at \(o\) must preserve both
eigenspaces. Thus
\[
 D_o\operatorname{Isom}(M,g)\subseteq O(V_p)\times O(V_q).
\]
The derivative map is injective. In fact, if \(\varphi(o)=o\) and
\(D_o\varphi=I\), then
\(\varphi(\exp_o v)=\exp_o(D_o\varphi v)=\exp_o v\); completeness and
Hopf--Rinow make \(\exp_o\) surjective. The given action already realizes
every element of \(O(V_p)\times O(V_q)\) as a derivative. Hence every
isometry has the same derivative as a unique element of the given group,
and injectivity identifies them. \(\square\)

### Corollary 2.2 (equal factor dimensions)

For \(k\ge2\) and \(0<\theta<1\),
\[
 \operatorname{Isom}(\mathcal S^{k,k}_\theta)=O(k)\times O(k)
\]
as the given group of transformations. In particular no additional factor
exchange isometry is possible.

**Proof.** Here
\[
 \rho_p-\rho_q=(k-1)(\lambda_p-\lambda_q)
              =(k-1)(\theta-1)\lambda_q\ne0.
\]
Positive curvature gives \(\lambda_q>0\). Apply Proposition 2.1.
\(\square\)

For unequal \(p,q\), the Ricci condition is instead
\((p-1)\lambda_p-(q-1)\lambda_q+(q-p)\mu\ne0\).
The prescribed ratio alone need not verify this condition. The full
curvature operator removes this limitation.

### Proposition 2.3 (curvature determination in all factor dimensions)

Under the geometric hypotheses of Proposition 2.1, with \(p,q\ge2\),
assume only that \(\lambda_p\ne\lambda_q\), rather than that
\(\rho_p\ne\rho_q\). Then the full isometry group is the given
\(O(p)\times O(q)\).

**Proof.** The first part of Proposition 2.1 shows that every isometry
fixes \(o\), independently of any eigenvalue assumption. Set
\(E=T_oM=V_p\oplus V_q\). Choose an orthonormal basis adapted to this
decomposition. Reflection of any single basis vector belongs to the given
isotropy group. For any two different unordered coordinate pairs, a
reflection in an index belonging to exactly one pair changes the sign of
one corresponding basis bivector and not the other. Invariance therefore
forces every off-diagonal coefficient of the curvature operator
\(\mathcal R:\Lambda^2E\to\Lambda^2E\) to vanish. Permutations within
each summand identify the diagonal entries, giving
\[
 \mathcal R=
 \lambda_p I_{\Lambda^2V_p}
 \oplus\mu I_{V_p\wedge V_q}
 \oplus\lambda_q I_{\Lambda^2V_q}.
\]
This argument uses the full orthogonal groups, including reflections,
and applies without exceptions when a factor has dimension two or four.

Because \(\lambda_p\ne\lambda_q\), at least one of these two pure
curvature values differs from \(\mu\). Its full eigenspace is thus
\(\Lambda^2V\) for one of the two factors \(V\). This eigenspace
determines \(V\) intrinsically. Indeed, identify bivectors with
skew-adjoint endomorphisms by
\[
 A_{u\wedge v}(x)=\langle u,x\rangle v-\langle v,x\rangle u.
\]
Then
\[
 \operatorname{span}\{A_\omega x:
       \omega\in\Lambda^2V,\ x\in E\}=V.
\]
The inclusion into \(V\) is immediate; the reverse inclusion follows by
choosing, for each nonzero \(v\in V\), a unit vector \(u\in V\) orthogonal
to \(v\), since \(A_{u\wedge v}u=v\). Here \(\dim V\ge2\) is essential.

The derivative of any isometry at \(o\) commutes with \(\mathcal R\),
so it preserves this eigenspace and its intrinsically defined support
\(V\). Orthogonality also preserves \(V^\perp\). Hence every such
derivative lies in \(O(V_p)\times O(V_q)\). The given group realizes
all of these derivatives, and the injectivity argument in Proposition
2.1 identifies every isometry with its corresponding group element.
\(\square\)

### Corollary 2.4 (the entire asymmetric Lavoyer--Peachey family)

For all \(p,q\ge2\) and \(0<\theta<1\),
\[
 \operatorname{Isom}(\mathcal S^{p,q}_\theta)=O(p)\times O(q).
\]

**Proof.** The prescribed tip curvatures satisfy
\(\lambda_p=\theta\lambda_q\) and \(\lambda_q>0\), so they are distinct.
Apply Proposition 2.3. \(\square\)

This extends Corollary 2.2 to unequal factor dimensions, even if the
tip Ricci eigenvalues happen to coincide. It does not establish the
strict total-curvature gap or remove the external input in Section 3.
No claim of literature priority is made.

## 3. A four-dimensional collapse consequence of recent literature

**External input [F, Theorem 1.1].** Fang's preprint, submitted 27 October
2026, states that a complete connected four-dimensional steady gradient
Ricci soliton with positive curvature operator satisfying (1) admits an
effective isometric \(O(3)\) action preserving its potential, with principal
orbits diffeomorphic to \(S^2\).

The statement and its noncollapse convention were checked in the primary
source, arXiv:2610.36563v1. The long analytic proof of that preprint is an
external dependency, not a proof supplied or independently re-established
here.

### Corollary 3.1 (conditional on the cited symmetry theorem)

Assuming [F, Theorem 1.1], no \(\mathcal S^{2,2}_\theta\), \(0<\theta<1\),
satisfies (1) for any \(\kappa>0\). Equivalently, there exist points \(x_j\)
and radii \(r_j>0\) with
\[
 \sup_{B_g(x_j,r_j)}|\mathrm{Rm}|\le r_j^{-2},
 \qquad
 \frac{\operatorname{Vol}_g B_g(x_j,r_j)}{r_j^4}\longrightarrow0.
\]

**Proof.** If (1) held, [F, Theorem 1.1] would provide an effective
\(O(3)\) action by isometries. Corollary 2.2 identifies the entire isometry
group with \(O(2)\times O(2)\). An effective \(SO(3)\) subgroup would therefore
give an injective Lie algebra map
\[
 \mathfrak{so}(3)\longrightarrow
 \mathfrak{so}(2)\oplus\mathfrak{so}(2),
\]
which is impossible: the source has dimension three, whereas the target
has dimension two (and is abelian). The asserted sequence follows by
negating (1) for \(\kappa=1/j\). \(\square\)

This rules out using all-scale \(\kappa\)-noncollapse as a property to be
proved for the four-dimensional family, subject to the explicit external
input. It does **not** rule out \(\alpha>0\). The curvature-controlled balls
in (1) and the distance blow-down defining the ideal boundary are different
geometric scales.

## 4. The actual reduction has two independent variables

The principal orbit is
\[
 (O(p)\times O(q))/(O(p-1)\times O(q-1))
 \cong S^{p-1}\times S^{q-1},
\]
of dimension \(p+q-2\). Thus the full soliton has **cohomogeneity two**,
including when \(p=q=2\). Its invariant hypersurfaces may have
cohomogeneity one; that does not turn the ambient equation into an ODE.

Choose one coordinate axis in each Euclidean summand and let
\[
 S=\operatorname{Fix}(O(p-1)\times O(q-1)).
\]
Here each subgroup fixes the chosen axis pointwise. In the standard model,
\(S\cong\mathbb R^2\). It is a closed, complete, totally geodesic surface.
At a principal point, its tangent plane is the fixed subspace of the
isotropy group; the orbit tangent representations have no fixed vectors.
This proves orthogonality to the principal orbit. The two sign changes
along the chosen axes act on \(S\) as commuting reflections. A quadrant
\(Q\) in \(S\) is the orbit space.

On the principal region, invariance under the separate orthogonal groups
therefore gives
\[
 g=h+a^2g_{S^m}+b^2g_{S^n},
 \qquad m=p-1,\quad n=q-1,
 \tag{4}
\]
where \(h\) is a metric on the two-dimensional quadrant, \(a,b>0\) there,
and \(a,b,F\) are functions of two variables. The potential can be taken
invariant: averaging preserves the steady equation; alternatively the
normalized potential is unique since a difference of two such potentials
has parallel gradient, which \(\operatorname{Ric}>0\) excludes.

All derivatives in the following equations are for \(h\):
\[
 \operatorname{Ric}_h-m\frac{\nabla^2a}{a}
                       -n\frac{\nabla^2b}{b}=\nabla^2F,
 \tag{5}
\]
\[
 \Delta a+\langle\nabla F,\nabla a\rangle
      +n\frac{\langle\nabla a,\nabla b\rangle}{b}
    =\frac{m-1}{a}(1-|\nabla a|^2),
 \tag{6}
\]
\[
 \Delta b+\langle\nabla F,\nabla b\rangle
      +m\frac{\langle\nabla a,\nabla b\rangle}{a}
    =\frac{n-1}{b}(1-|\nabla b|^2).
 \tag{7}
\]

For completeness, if \(X\) is a unit base vector and \(U,V\) are unit
vectors in the two fibers, the connection formula
\(\nabla_U X=X(\log a)U\) and its \(b\)-analogue give
\[
 K(X,U)=-\frac{\nabla^2a(X,X)}a,
 \quad K(X,V)=-\frac{\nabla^2b(X,X)}b,
 \quad K(U,V)=-\frac{\langle\nabla a,\nabla b\rangle}{ab}.
 \tag{8}
\]
Tracing the first identities gives (5). The sectional curvatures within
the first fiber are \((1-|\nabla a|^2)/a^2\); tracing these and (8), then
equating to \(\nabla^2F(U,U)=\langle\nabla F,\nabla a\rangle/a\), gives
(6). The second fiber gives (7). When a fiber is one-dimensional there
are no within-fiber planes, explaining the vanishing coefficient.

### Proposition 4.1 (exact weighted identities)

Put \(V=a^mb^n\) and \(\Delta_F=\Delta+\langle\nabla F,\nabla\cdot\rangle\).
Every solution of (6)--(7) satisfies
\[
 \Delta_F V=
 V\left(\frac{m(m-1)}{a^2}+\frac{n(n-1)}{b^2}\right),
 \tag{9}
\]
\[
 \operatorname{div}\left(e^F V\nabla\log\frac ab\right)
 =e^F V\left(\frac{m-1}{a^2}-\frac{n-1}{b^2}\right).
 \tag{10}
\]
In particular, for \(p=q=2\),
\[
 \operatorname{div}(e^F\nabla(ab))=0,
 \qquad
 \operatorname{div}\left(e^Fab\nabla\log\frac ab\right)=0.
 \tag{11}
\]

**Proof.** Let \(u=\log a\), \(v=\log b\), and
\(\mathcal L=\Delta_F+\langle\nabla\log V,\nabla\cdot\rangle\).
Dividing (6) by \(a\), and (7) by \(b\), gives
\[
 \mathcal Lu=(m-1)a^{-2},\qquad
 \mathcal Lv=(n-1)b^{-2}.
\]
Subtracting proves (10). Multiplying the first identity by \(m\) and the
second by \(n\), and using
\(V^{-1}\Delta_FV=\Delta_F\log V+|\nabla\log V|^2\), proves (9).
Equation (11) follows by putting \(m=n=1\). \(\square\)

These are identities on the open quadrant. The weights and logarithms
degenerate at its axes, so an integration argument requires justified
boundary terms. They alone provide no estimate at infinity for \(\alpha\).

## 5. A flat-base reduction forces a product

### Proposition 5.1 (flat-section rigidity)

Let \((M^{p+q},g)\) be complete with the standard \(O(p)\times O(q)\)
symmetry and \(\sec_g\ge0\). If the metric \(h\) in (4) is flat, then
\((M,g)\) is a product of complete rotational metrics on
\(\mathbb R^p\) and \(\mathbb R^q\).

**Proof.** The complete fixed surface \(S\cong\mathbb R^2\) above is flat
on the open quadrants and hence, by smoothness, everywhere. Being complete
and simply connected, it is the Euclidean plane. The Weyl reflections
are Euclidean reflections whose axes meet orthogonally at \(o\). Thus the
quotient is the whole Euclidean quadrant \(x,y\ge0\), and
\(h=dx^2+dy^2\).

Fix \(x>0\). Smooth reflection across \(y=0\) gives
\(a_y(x,0)=0\). Equation (8) and \(\sec_g\ge0\) give
\(a_{yy}(x,y)\le0\) for \(y>0\). Hence \(a_y(x,y)\le0\). If it were
negative anywhere, concavity would keep it bounded above by a negative
constant thereafter, forcing \(a(x,y)<0\) for sufficiently large \(y\).
But \(a(x,y)>0\) for all \(y\ge0\), by the standard orbit types. Thus
\(a_y\equiv0\). The same argument gives \(b_x\equiv0\), so
\[
 g=(dx^2+a(x)^2g_{S^{p-1}})
   +(dy^2+b(y)^2g_{S^{q-1}}).
\]
Smoothness at the axes closes each factor to Euclidean topology;
completeness is inherited from the product. \(\square\)

### Corollary 5.2 (the four-dimensional steady case)

If additionally \(p=q=2\), \(\operatorname{Ric}>0\), and
\(\operatorname{Ric}=\nabla^2F\), the two factors in Proposition 5.1 are
scaled Hamilton cigars.

**Proof.** Use an invariant potential. In the product coordinates the
mixed Hessian equation gives \(F_{xy}=0\), hence \(F=F_1(x)+F_2(y)\).
For a factor \(ds^2+w(s)^2d\phi^2\), smoothness gives
\(w(0)=0,w'(0)=1\). The two steady equations are
\[
 -\frac{w''}{w}=F_1'',\qquad -w''=F_1'w'.
\]
Combining them gives \((F_1'/w)'=0\), so \(F_1'=c w\), and
\[
 w'=1-\frac c2w^2.
\]
Positive Gaussian curvature \(-w''/w=cw'>0\) near the origin gives
\(c>0\). With \(A=\sqrt{2/c}\), the solution is
\[
 w(s)=A\tanh(s/A),\qquad
 F_1(s)=2\log\cosh(s/A)+\mathrm{constant}.
\]
This is precisely a scaled cigar. Apply the same calculation to the other
factor. \(\square\)

### Exact counterexample after weakening strict curvature

For any \(A,B>0\), the metric and potential
\[
 g=ds^2+A^2\tanh^2(s/A)d\phi^2
   +dt^2+B^2\tanh^2(t/B)d\psi^2,
\]
\[
 F=2\log\cosh(s/A)+2\log\cosh(t/B)
\]
give a complete standard \(O(2)\times O(2)\)-symmetric steady soliton with
positive Ricci tensor and nonnegative curvature operator. At the vertex,
\[
 \lambda_p=2/A^2,\quad\lambda_q=2/B^2,
 \quad\lambda_p/\lambda_q=B^2/A^2.
\]
Thus any prescribed \(0<\theta<1\) is obtained by \(B=\sqrt\theta A\).
Nevertheless both compact-group orbits collapse at distance scale: the
two angular radii are bounded. The projection to the Euclidean quadrant
is distance nonincreasing, and choosing angular representatives changes
distances by at most \(\pi A+\pi B\). After division by the distance scale
this error vanishes. The asymptotic cone is therefore the Euclidean
quadrant, its ideal boundary is \([0,\pi/2]\), and \(\alpha=\beta=0\).

The example is not a counterexample to the question in [LP]: every mixed
sectional curvature is zero, so its curvature operator is not positive.
It proves that positivity of Ricci, nonnegativity of the curvature
operator, the tip-curvature ratio, completeness, and symmetry do not
suffice to prove \(\alpha>0\). Proposition 5.1 also explains why replacing
the two-variable base metric by a flat quadrant cannot capture the
strictly positively curved target family.

## 6. An intrinsic criterion for boundary-orbit noncollapse

Let \(A:M\to[0,\infty)\) be the radius of the \(O(p)\)-orbit: on the
principal region it is the function \(a\) in (4). It vanishes precisely on
\(\Sigma_q=\operatorname{Fix}(O(p))\). Write \(u\in S^{p-1}\) for the
direction in the first Euclidean summand of the standard action.

### Lemma 6.1 (equivariant Lipschitz radius map)

On any complete standard \(O(p)\times O(q)\)-symmetric manifold with
\(\sec\ge0\), the map
\[
 \Phi_p:M\longrightarrow\mathbb R^p,
 \qquad \Phi_p=A u\quad\hbox{off }\Sigma_q,
 \qquad \Phi_p=0\quad\hbox{on }\Sigma_q,
 \tag{12}
\]
is 1-Lipschitz.

**Proof.** Fix the subgroup \(O(p-1)\) fixing a chosen axis in the first
summand. Its fixed set \(N\cong\mathbb R\times\mathbb R^q\) is complete
and totally geodesic. Its two halves are interchanged by sign change in
the first coordinate, with common fixed hypersurface \(\Sigma_q\).
On the positive half \(N_+\), the ambient metric has the warped-product
form
\[
 g=g_N+A^2g_{S^{p-1}}.
\]
Equation (8), now with the whole \(N_+\) as base, shows
\(\nabla_N^2 A\le0\). Smoothness of the standard spherical collapse
extends \(A\) to a signed smooth function across \(\Sigma_q\) in \(N\),
with \(|\nabla_N A|=1\) there.

Suppose \(|\nabla_N A|(z)=k>1\) at a point in \(N_+\). Take the complete
unit-speed geodesic with initial velocity \(-\nabla_N A/k\). While it stays
in \(N_+\), concavity gives \((A\circ\gamma)'\le-k\). It must therefore
first meet \(\Sigma_q\) at a time \(T\le A(z)/k\). At that first meeting,
smoothness and \(|\nabla_N A|=1\) imply
\((A\circ\gamma)'(T)\ge-1\), contradicting its limit from the left,
which is at most \(-k\). Hence \(|\nabla_N A|\le1\).

For a tangent vector decomposed into its \(N\) component and spherical
component, the squared norm of its image under (12) is
\((dA)^2+A^2|du|^2\), whereas its squared ambient norm is
\(|v_N|^2+A^2|du|^2\). The gradient bound proves the differential estimate.
Smooth standard collapse makes (12) extend smoothly across \(\Sigma_q\)
(equivalently, the signed radius is odd with the required first-order
normal coefficient). Integrating along curves proves the global bound.
\(\square\)

The restriction to the other fixed submanifold has the intrinsic form
\[
 (\Sigma_p,g|_{\Sigma_p})
   =([0,\infty)\times S^{p-1},ds^2+w(s)^2g_{S^{p-1}}),
 \tag{13}
\]
with the sphere collapsed at zero. Smoothness and nonnegative radial
curvature give
\(w(0)=0,w'(0)=1,w''\le0\). Completeness makes the radial interval
infinite; the standard orbit types give \(w>0\) away from zero. A negative
derivative would force \(w\) eventually negative, so \(w'\ge0\). Thus
\[
 c_p:=\lim_{s\to\infty}w'(s)
      =\lim_{s\to\infty}\frac{w(s)}s\in[0,1].
 \tag{14}
\]

### Proposition 6.2 (collapse criterion on the fixed submanifold)

Let \((M,g,F)\) be a complete connected steady soliton with the standard
\(O(p)\times O(q)\) symmetry, \(\sec\ge0\), \(\operatorname{Ric}>0\), and
a critical point at the common fixed point \(o\). The \(O(p)\)-action on
the asymptotic cone is nontrivial if and only if \(c_p>0\).
For the family in [LP], this says exactly
\[
 \boxed{\alpha>0\quad\Longleftrightarrow\quad c_p>0.}
 \tag{15}
\]

**Proof.** First suppose \(c_p>0\). Take opposite points \(y_s,y_s^-\) in
the \(O(p)\)-orbit of intrinsic radius \(s\) in (13). They lie in
\(B_g(o,s)\), since their intrinsic radial paths have length \(s\).
Lemma 6.1 gives
\[
 d_g(y_s,y_s^-)\ge|\Phi_p(y_s)-\Phi_p(y_s^-)|
                 =2w(s)\ge2c_p s.
 \tag{16}
\]
After distance rescaling by \(s^{-1}\), the isometry \(-I_p\) therefore
has displacement at least \(2c_p\) somewhere in the unit ball. Compactness
of bounded balls in the blow-down limit, with the induced compact-group
action, shows this action cannot be trivial. For a nonnegatively curved
manifold these blow-downs are the asymptotic cone, with the natural action
induced by isometries fixing \(o\).

For the converse suppose \(c_p=0\). We prove the uniform estimate
\[
 \sup_{x\in B_g(o,r)}A(x)=o(r).
 \tag{17}
\]
For \(x\notin\Sigma_q\), choose a point \(y\in\Sigma_p\) nearest to
\(x\). Such a point exists because \(\Sigma_p\) is closed and \(M\) is
complete. If \(x\in\Sigma_p\) the next comparison is immediate. Otherwise
take a minimizing geodesic \(\gamma\) from \(y\) to \(x\), perpendicular
to \(\Sigma_p\) at its initial point. Here \(y\ne o\): if \(y=o\), its
initial vector would lie in \(T_o\Sigma_q\), forcing the whole geodesic
into the totally geodesic \(\Sigma_q\).

The subgroup \(O(p-1)\) fixing the first-summand direction of \(y\)
fixes both \(y\) and every normal vector to \(\Sigma_p\) there, so
\(\gamma\) lies in the corresponding complete \(N\). It cannot cross
\(\Sigma_q\). Indeed, at a first meeting, reflecting the initial segment
by \(-I_p\) leaves the meeting point fixed and moves its initial point
within \(\Sigma_p\). It gives another path of the minimum possible length
from \(\Sigma_p\) to \(x\), with a corner unless \(\gamma'\) is tangent to
\(\Sigma_q\) at that meeting. The corner contradicts length minimality;
tangency contradicts geodesic uniqueness and total geodesicity, since
\(y\notin\Sigma_q\). Thus \(\gamma\) stays in \(N_+\).

Along this geodesic, \(A''\le0\). At \(y\), \(A'(0)=0\): the normal
space to \(\Sigma_p\) is the standard \(O(q)\) representation, and the
differential of the invariant function \(A\) vanishes on it. Consequently
\[
 A(x)\le A(y)=w(s_y),
 \tag{18}
\]
where \(s_y\) is the intrinsic radial coordinate in \(\Sigma_p\).

It remains to compare this intrinsic coordinate with ambient distance;
no global convexity of \(\Sigma_p\) is assumed. Normalize \(F(o)=0\).
By (3), \(|\nabla F|\le H:=\sqrt{R(o)}\). On an intrinsic radial geodesic
of \(\Sigma_p\),
\[
 F''=\operatorname{Ric}(\partial_s,\partial_s)>0,
 \qquad F'(0)=0.
\]
Thus there are constants \(\varepsilon>0,C_0\) with
\(F(y)\ge\varepsilon s_y-C_0\), uniformly in direction by symmetry.
On the other hand \(F(y)\le H d_g(o,y)\). Nearest-point minimality gives
\(d_g(x,y)\le d_g(x,o)\), hence \(d_g(o,y)\le2d_g(o,x)\). Therefore
\[
 s_y\le C_1d_g(o,x)+C_2.
\]
Together with (18), monotonicity of \(w\), and \(w(s)=o(s)\), this proves
(17). Points of \(\Sigma_q\) have \(A=0\) and cause no exception.

The ambient diameter of an \(O(p)\)-orbit is at most its intrinsic
diameter \(\pi A\). Estimate (17) therefore makes the whole \(O(p)\)
action converge to the identity uniformly on each bounded blow-down
ball. Its action on the asymptotic cone is trivial. Finally, in [LP]'s
warped boundary description the first group acts trivially precisely
when \(a\equiv0\), equivalently \(\alpha=0\). This proves (15).
\(\square\)

### Corollary 6.3 (the two-dimensional fixed-surface test)

Under the hypotheses of Proposition 6.2, if \(p=2\), then
\[
 c_p=1-\frac1{2\pi}\int_{\Sigma_p}K_{\Sigma_p}\,dA.
\]
In particular, for \(\mathcal S^{2,q}_\theta\),
\[
 \boxed{\alpha>0\quad\Longleftrightarrow\quad
 \int_{\Sigma_p}K_{\Sigma_p}\,dA<2\pi.}
 \tag{19}
\]

**Proof.** In (13), \(K=-w''/w\) and \(dA=w\,ds\,d\phi\). Thus
\[
 \int_{\Sigma_p}K\,dA
 =-2\pi\int_0^\infty w''(s)\,ds=2\pi(1-c_p).
\]
Apply Proposition 6.2. \(\square\)

### Proposition 6.4 (actual-axis integral and curvature alternatives)

Let \(p=q=2\) and suppose the complete soliton in Proposition 6.2
has positive curvature operator. On \(\Sigma_p\), write
\(f=F_0'\), \(k=-w''/w\), and let \(U,V,N\) be the sectional
curvatures of a radial--normal, angular--normal, and normal
two-plane, respectively. Thus \(U\) and \(V\) are generally
different. The axis equations and scalar identity are
\[
 f'=k+2U,\qquad f\,w'/w=k+2V,\qquad
 C=f^2+2k+4U+4V+2N,\quad C=R(o).
\]
There is a limit \(f\to\ell\in(0,\sqrt C]\), and
\[
 \int_0^\infty(k+2U)\,ds=\ell,\qquad
 w\text{ is bounded}\ \Longleftrightarrow\
 \int_1^\infty V\,ds<\infty.
\]
If \(w\to\infty\), then
\[
 \frac{2\int_1^s V(r)\,dr}{\log w(s)}\longrightarrow\ell.
\]
Furthermore \(k,U,V\to0\) and \(N\to(C-\ell^2)/2\).

**Proof.** Strict positivity gives \(f'>0\), while the Hamilton
identity bounds \(f\le\sqrt C\). Since \(f(0)=0\), its limit is
positive and integration gives the first assertion.
Completeness, positivity of \(w\), strict concavity, and
\(w(0)=0\) imply \(0<w'/w\le1/s\).
Consequently \(0<k+2V\le\sqrt C/s\), so \(k,V\to0\).
Integrating the angular equation from \(1\) to \(s\) gives
\[
 2\int_1^s V\,dr
 =\int_1^s f(r)\,d\log w(r)-\int_1^s k(r)\,dr.
\]
The last integral is bounded, and \(f(1)\le f(r)\le\ell\).
This proves the boundedness equivalence. If \(\log w(s)\to\infty\),
the weighted average of \(f(r)\to\ell\) proves the stated limit.

For the remaining pointwise limits, the standard curvature
derivative estimate is needed, not just integrability.
Positive curvature operator and \(R\le C\) bound
\(|\mathrm{Rm}|\) globally. The complete steady flow is eternal
and has the same curvature bound at every time. Applying Shi's
gradient estimate [S] on a fixed past time interval bounds
\(|\nabla\mathrm{Rm}|\) globally. Along the unit radial geodesic,
\(f'=\operatorname{Ric}(\partial_s,\partial_s)\) is therefore
uniformly Lipschitz. A nonnegative integrable uniformly Lipschitz
function tends to zero: otherwise disjoint intervals of a fixed
positive length would carry a fixed positive integral.
Hence \(f'\to0\), so \(U\to0\). The scalar identity now gives
the limit of \(N\). \(\square\)

In particular, an unbounded sublinear axis has nonsummable \(V\)
but summable \(U\). It cannot satisfy \(U=V\), or even an eventual
bound \(U\ge\varepsilon V\) with fixed \(\varepsilon>0\).
These statements apply to actual complete axes; they do not prove
\(c_p>0\). If \(C=\ell^2\), the limiting normal curvature vanishes,
so even a nonflat transverse cigar limit has not been supplied.
If \(C>\ell^2\), the curvature limit alone supplies neither the
product-end expansion nor the mode coefficients used in Section 7.

### Example 6.5 (unbounded sublinear data pass the axis tests)

The preceding necessary identities do not force a bounded radius,
even with every axis sectional curvature positive and a prescribed
tip ratio. Fix \(0<\theta<1\), choose \(c>6\), and take
\[
 0<a<\min\left\{\frac12,\,
          \frac{48(\theta^{-1}-1)}{c^2\pi^2}\right\},\qquad B=ca.
\]
Define
\[
 w(s)=\frac{s}{(1+s^2)^a},\quad
 F_0(s)=B\left(s\arctan s-\tfrac12\log(1+s^2)\right),\quad
 f=B\arctan s,
\]
and
\[
\begin{split}
 k&=\frac{2a(3+(1-2a)s^2)}{(1+s^2)^2},\\
 U&=(f'-k)/2,\qquad V=(fw'/w-k)/2,\\
 C&=4B+12a(\theta^{-1}-1),\\
 N&=\tfrac12(C-f^2-2f'-2fw'/w+2k).
\end{split}
\]
All expressions have their smooth axis limits at \(s=0\);
\(w\) is odd analytic and \(F_0,k,U,V,N\) are even analytic.
Directly,
\[
 w'=\frac{1+(1-2a)s^2}{(1+s^2)^{a+1}}>0,\qquad k=-w''/w>0.
\]
The estimates
\[
 k\le\frac{6a}{1+s^2},\qquad
 \frac{k}{w'/w}
 =\frac{2as(3+(1-2a)s^2)}
        {(1+s^2)(1+(1-2a)s^2)}
 \le\frac{6as}{1+s^2}\le6a\arctan s
\]
prove \(U,V>0\), since \(B>6a\).
Also \(f'\le B\), \(fw'/w\le B\), and \(f^2\le B^2\pi^2/4\),
so
\[
 N\ge6a(\theta^{-1}-1)-B^2\pi^2/8>0.
\]
At the common vertex,
\[
 k(0)=6a,\quad U(0)=V(0)=(B-6a)/2,\quad N(0)=6a/\theta.
\]
Thus the required pure-curvature ratio is \(\theta\).
The defining formulas verify all three axis identities above.
Nevertheless \(w(s)\sim s^{1-2a}\to\infty\) and \(w'(s)\to0\).

These are compatible analytic axis data, not a complete soliton
and not a counterexample to the LP conjecture. No full ambient
existence or continuation assertion is made here. They show
precisely why the axis identities and the exclusion of the special
bounded cigar profile do not eliminate the unbounded sublinear
branch of \(c_p=0\).


### Proposition 6.6 (ambient limits along every unbounded axis)

Under the hypotheses of Proposition 6.4, assume \(w(s)\to\infty\),
without assuming \(w(s)/s\) has a positive limit. Set
\(c=(C-\ell^2)/2\). For every sequence \(p_i\in\Sigma_p\) with
\(s(p_i)\to\infty\), the unrescaled pointed manifolds have smooth
subsequential limits. Every such limit, with potential normalized
to vanish at its basepoint, is
\[
 \begin{cases}
 \mathbb R^2\times\mathrm{Cigar}_c,\quad
 F_\infty=\ell x+2\log\cosh(\sqrt{c/2}\,t),&c>0,\\
 \mathbb R^4,\quad F_\infty=\ell x,&c=0.
 \end{cases}
\]
The cigar has metric
\(dt^2+(2/c)\tanh^2(\sqrt{c/2}\,t)d\psi^2\), and the basepoint
is at its tip. These are actual products, not merely statements
about universal covers.

**Proof.** The Hamilton identity bounds the curvature tensor and
\(|\nabla F|\). The steady flow is complete and eternal; Shi's
estimates give uniform bounds on all curvature derivatives.
The Gromoll--Meyer injectivity estimate, in the explicit formulation
[CZ, Theorem 4.2.1], gives
\(\operatorname{inj}(M,g)\ge\pi/\sqrt{\sup\sec_g}>0\).
Smooth pointed compactness, equivalently [CZ, Theorem 4.1.5],
therefore gives complete subsequential limits.
The functions \(F-F(p_i)\) converge as well, since their gradients
are bounded and their Hessians are Ricci. The limit is a complete
steady soliton with nonnegative curvature operator.
At its basepoint Proposition 6.4 gives
\[
 \operatorname{Ric}_\infty=\operatorname{diag}(0,0,c,c),
 \qquad \nabla F_\infty=\ell e_s.
\]
The normal \(SO(2)\) isotropy also converges: it fixes the basepoints
and retains the identity representation on the first two tangent
directions and the standard rotations on the normal plane.

We first rule out a compact quotient in the angular flat direction.
Let \(\mathcal A\) be the first orbit radius and \(K=\partial_\phi\)
the period-\(2\pi\) Killing field. Lemma 6.1 and the warped metric give
\[
 |K|=\mathcal A,\qquad
 |\nabla K|^2=2|\nabla\mathcal A|^2\le2.
\]
With \(w_i=\mathcal A(p_i)\to\infty\), the fields \(X_i=K/w_i\)
converge on fixed balls to a unit parallel field \(X_\infty\).
Indeed \(|\mathcal A-w_i|\le r\) on \(B(p_i,r)\) and
\(|\nabla X_i|\le\sqrt2/w_i\); higher convergence follows from
the Killing-field identities and curvature derivative bounds.

This parallel field has an exact dual. The first summand direction
defines \(u:M\setminus\Sigma_q\to S^1\).
Along a minimizing path of length at most \(r\) from \(p_i\),
\(\mathcal A\ge w_i-r\) and
\(\int|d\phi|\le r/(w_i-r)\).
Hence the entire ball has its angular image in a single proper arc.
The real branches
\[
 y_i=w_i\arg(u/u(p_i)),\qquad
 dy_i=(w_i^2/\mathcal A^2)X_i^\flat
\]
agree on nested fixed balls. Their limit is one globally defined
function \(y\) with \(dy=X_\infty^\flat\),
\(|\nabla y|=1\), and \(\nabla^2y=0\).
Completeness gives a global \(\mathbb R_y\) product factor.
Also \(\langle\nabla y,\nabla F_\infty\rangle=0\).

If \(c>0\), the tensor strong maximum principle for the complete
bounded-curvature eternal limit [H; CZ, Section 2.2] makes the
Ricci nullspace parallel, of rank two.
The universal cover splits as \(\mathbb R^2\times N^2\).
The lifted potential is \(\ell x+G\), with \(\nabla G=0\)
at the surface basepoint. The surviving normal rotations act on
\(N^2\). Its polar steady equations are
\(G''=G'B'/B=-B''/B\).
Smoothness and the tip curvature \(c\) give \(G'=cB\) and
\(B'=1-cB^2/2\), yielding precisely the displayed complete cigar.
The polar identification extends globally on the simply connected
complete surface.

A deck transformation preserves \(y\), so its Euclidean part is
\((x,y)\mapsto(\pm x+a,y)\).
Every cigar isometry fixes its unique curvature maximum and
preserves \(G\). Preservation of the global potential
\(\ell x+G\), with \(\ell>0\), forces the plus sign and \(a=0\).
The transformation would then fix points over the cigar tip.
Deck transformations act freely, so the deck group is trivial.

If \(c=0\), the nonnegative scalar curvature vanishes at the
basepoint. The equation
\[
 (\Delta+\nabla F_\infty\cdot\nabla)R_\infty
       =-2|\operatorname{Ric}_\infty|^2
\]
and the strong maximum principle make the limit flat.
Now \(x=F_\infty/\ell\) and \(y\) are global orthogonal affine
unit coordinates, giving \(\mathbb R^2\times L^2\) with \(L^2\)
complete and flat. Normal rotations preserve both coordinates,
fix the basepoint, and act on \(L^2\) by all plane rotations.
Lift these rotations to the Euclidean universal cover fixing the
origin. They normalize the discrete deck group, and connectedness
makes their conjugation action on each deck transformation constant.
Every deck transformation thus commutes with all origin rotations,
and consequently fixes the origin. Freeness again makes the deck
group trivial. Thus \(L^2=\mathbb R^2\). \(\square\)

This is fixed-radius ambient convergence, including for unbounded
sublinear axes. It gives neither a convergence rate nor control
on transverse distances increasing with \(s\).
In particular it does not exclude \(w(s)/s\to0\).
When \(c=0\), it says nothing about limits after rescaling by a
curvature scale that tends to infinity.

### An obstruction to obtaining relative derivative bounds from axis tests

The axis identities in Proposition 6.4, even together with bounded
derivatives of every fixed order of the axis curvature functions, do
not imply the scale-dependent bounds
`|a^(j)|=O(a^(j+1))`, where `a=w'/w`. This matters for replacing the
expansion hypothesis in Proposition 6.7 by a general slowly varying
normal development.

Here is an explicit end counterexample. Fix `0<mu<1`, `ell,c>0`,
`epsilon!=0`, and `b>mu(1-mu)`. On a sufficiently late half-line put

    delta=epsilon*s^(-5/2)*cos(s),
    w=s^mu*(1+delta),  f=ell-b/s,  C=ell²+2c,
    k=-w''/w,  U=(f'-k)/2,  V=(f*w'/w-k)/2,
    N=(C-f²-2f'-2f*w'/w+2k)/2.

Then `w>0`, `w'>0`, `w''<0`, `w->infinity`, and `w'->0`. The three
axis identities hold exactly by substitution. Since every fixed
derivative of delta is `O(s^(-5/2))`, direct differentiation gives

    k=mu(1-mu)/s²+O(s^(-5/2)),
    U=(b-mu(1-mu))/(2s²)+O(s^(-5/2)),
    V=ell*mu/(2s)+O(s^(-2)),
    N=c+ell*(b-mu)/s+O(s^(-2)).

All four functions are strictly positive eventually, and every fixed
derivative is bounded. For example,
`w''/w=mu(mu-1)/s²+(delta''+2mu*delta'/s)/(1+delta)` proves the
first estimate with a strictly smaller error than its leading term.
The other estimates follow from the displayed exact definitions.

Nevertheless

    a''/a³ = 2/mu²+(epsilon/mu³)*s^(1/2)*sin(s)+O(s^(-1/2)).

Indeed `a=mu/s+delta'/(1+delta)` and
`a''=2mu/s³+epsilon*s^(-5/2)*sin(s)+O(s^(-7/2))`.
The ratio is unbounded with both signs along the two sequences where
sin(s)=1 and sin(s)=-1. Thus already `|a''|=O(a³)` fails.
These are end data satisfying the necessary axis tests, not an actual
ambient soliton. Bounded derivatives here concern only those axis
functions; no ambient Shi estimate is asserted for a constructed metric.
The example ends the proposed deduction from the established axis tests
and absolute derivative bounds. It does not prevent the full soliton
equations from imposing additional relative control.

### A conditional Riccati criterion for the strict slope gap

The ODE step of one possible general-profile argument has a precise
sufficient condition. Let w be a positive increasing concave C³ function
on a half-line, with `w->infinity`. Put

    a=w'/w,  u=log(w),  b=-a'/a².

Concavity gives `b>=1`. If, for a nonnegative integrable function e(u),

    db/du >= b²-1-e(u),

then `lim w'>0`. To prove this, set `z=b-1>=0` and
`E(u)=integral_u^infinity e(t)dt`. If z>E somewhere, then
`v=z-E` satisfies `v'>=2z+z²>=v²` while v>0, forcing finite-u blowup.
This contradicts smooth positivity of a at every finite u. Hence
`z<=E->0`. Integrating `z'-2z>=-e` to infinity now gives

    z(u)<=integral_u^infinity exp(-2(t-u))*e(t)dt,
    integral_(u0)^infinity z <= (1/2)*integral_(u0)^infinity e < infinity.

Finally `(log w')_u=1-b=-z`, so w' has a positive limit.
The change of variables also gives the exact identity

    a⁴+(a')²-a*a'' = a⁴*(1-b²+db/du).

Consequently a lower bound for this expression with an integrable
relative error would suffice for the strict gap. Such a bound has not
been derived from the actual LP normal curvature. The preceding
counterexample explains why absolute axis derivative estimates alone
cannot supply the relative expansion needed for that step.

### Proposition 6.7 (exclusion of analytic sublinear power-law axes)

An actual complete axis under the hypotheses of Proposition 6.4
cannot simultaneously have \(c>0\) and the expansion
\[
 w(s)=w_0s^\mu W(1/s),\qquad
 F_0'(s)=\ell-\frac b s+\frac1{s^2}G(1/s),
 \qquad 0<\mu<1,
\]
where \(w_0>0\), and \(W,G\) are real analytic near zero with
\(W(0)=1\). The obstruction is a negative sectional curvature
at normal distances \(t=\gamma\log s\), for a fixed sufficiently
large \(\gamma\). The proof first forces two exact coefficient
equalities, then shows that the meridional sectional curvature is
\[
 -\frac{\mu^2(1-\mu^2)(d^2+\ell^2)}{d^4}s^{-4}
       +o(s^{-4})<0,\qquad d=\sqrt{2c}.
\]
This excludes the whole analytic power-law sublinear regime with a
nonflat transverse limit. It does not assert that every sublinear
axis has such an expansion, or resolve the branch \(c=0\).

**Proof.** Write \(\varepsilon=1/s\), \(d=\sqrt{2c}>0\),
\(\omega=W'(0)\), \(g=G(0)\), and
\(N_1=\ell(b-\mu)\). The axis identities give
\[
 k=\mu(1-\mu)\varepsilon^2+O(\varepsilon^3),\quad
 U=\tfrac12(b-\mu(1-\mu))\varepsilon^2+O(\varepsilon^3),\quad
 V=\tfrac12\ell\mu\varepsilon+O(\varepsilon^2),\quad
 N=c+N_1\varepsilon+O(\varepsilon^2).
\]
Positivity of \(U\) proves the lower bound.

For the upper bound, normalize the normal development by
\[
 u=(\log H,\log(a/w),\log(b_{\rm orb}/B),F-F_0-\mathcal U),
\]
where \(B=(2/d)\tanh(dt/2)\), \(\mathcal U=2\log\cosh(dt/2)\),
and \(b_{\rm orb}\) denotes the second orbit radius.
The normalized logarithmic equations have analytic coefficients
in \(\varepsilon\), since
\(\partial_s=-\varepsilon^2\partial_\varepsilon\);
the nonanalytic values of \(\log w\) and \(F_0\) enter only
through their analytic derivatives. The general axis jets are
\[
 H=1-\tfrac12Ut^2+\cdots,\quad
 a/w=1-\tfrac12Vt^2+\cdots,\quad
 b_{\rm orb}/t=1-\tfrac16Nt^2+\cdots,\quad
 F-F_0=\tfrac12(U+V+N)t^2+\cdots.
\]
The parameter-dependent Fuchsian construction used in Section 7
therefore gives a common initial interval with
\(u=\varepsilon u_1+\varepsilon^2u_2+O(\varepsilon^3)\).
The actual steady soliton is real analytic [K, Corollary 1.3],
as is its potential in geodesic normal coordinates. Its fixed axis
and normal exponential gauge are analytic as well. Uniqueness of
the normal recurrence therefore identifies its germ with this
parameter-dependent solution.
In particular it does not replace \(V\) by \(U\).

Put \(D=d\coth(dt)\), \(L=d\operatorname{csch}(dt)\).
At first order \(h_1=0\), and the angular and second-radius
corrections satisfy
\[
 \alpha_1''+D\alpha_1'=-\ell\mu,\qquad
 \beta_1=N_1\left(-\frac1{d^2}+
                    \frac{t}{d\sinh(dt)}\right).
\]
Thus \(\alpha_1'\to-\ell\mu/d\),
\(\beta_1\to-N_1/d^2\), and \(\beta_1'\) decays exponentially
up to a polynomial factor. If \(Q_j\) is the sum of the four
components of \(u_j\), the second-order equations include
\[
 \beta_2''+D\beta_2'+LQ_2'
       =\ell\beta_1-\beta_1'Q_1',\qquad
 Q_2''=-2L\beta_2'-(\alpha_1')^2-(\beta_1')^2.
\]
The homogeneous normal subsystem is stable plus a constant mode.
These equations give polynomial growth of the second-order
coefficients, \(Q_2'=O(1+t)\), and
\[
 \beta_2'\longrightarrow-\ell N_1/d^3.
\]
Unknown second-order axis coefficients affect transient terms,
not this limit. For completeness, with \(\omega=W'(0)\), the other
second-order equations are
\[
 h_2''+Dh_2'=-[b-\mu(1-\mu)],\qquad
 \alpha_2''+D\alpha_2'
 =\mu(1-\mu+b)+\ell(\omega+\alpha_1)-\alpha_1'Q_1'.
\]
They give \(h_2=O(1+t)\), \(\alpha_2=O((1+t)^2)\),
\(\beta_2=O(1+t)\), and \(Q_2=O((1+t)^2)\);
hence every component of \(u_2\), including \(f_2\), is polynomially
bounded.

Here the expansion must again be controlled on a growing interval.
Fix \(\varepsilon_0>0\) and use \(z=\varepsilon/\varepsilon_0\)
on a fixed complex disc around \(1\).
Retain \(X=(h,u_t,Z)\), with \(Z=u_s\).
The spatial principal coefficients are \(O(\varepsilon_0)\),
because \(\partial_s=-\varepsilon_0z^2\partial_z\).
The zero-parameter linear system is
\[
\begin{split}
 p_h'&=-Dp_h+\ell Z_h,\qquad
 p_a'=-Dp_a-\ell Z_a,\\
 p_b'&=-Dp_b-Lq-\ell Z_b,\qquad q'=-2Lp_b,\\
 h'&=p_h,\qquad Z'=0,
\end{split}
\]
where \(q=p_h+p_a+p_b+p_f\).
Its \((p,Z)\) subsystem has a bounded forward propagator;
the \(h\) component has at most linear growth.
Equivalently its full propagator is bounded by \(C(1+t-\sigma)\).
This follows by solving the stable equations and using the
integrability of \(L,D-d\); the constant \(Z_b\) forcing enters
the stable \(p_b\) equation, while the \(q\) coupling is integrable.

On \(t\le T=\gamma\log(1/\varepsilon_0)\), use
\(C_T=C(1+T)\) as the propagator bound in Duhamel's analytic
majorant. The scalar transport argument of Proposition 7.11,
with its spatial and state majorants retained, now has source
\(O(C_T\varepsilon_0)\), derivative coefficient
\(O(C_T\varepsilon_0)\), and quadratic coefficient \(O(C_T)\).
Its amplitude is \(O(\varepsilon_0(1+T)^2)\);
its radius loss is \(O(\varepsilon_0(1+T)^2)\).
The nonlinear amplification also tends to zero, since
\(\varepsilon_0\) times any fixed polynomial in \(T\) tends to
zero.

For the second-order approximation retain its exact derivative
\[
 Z_{\rm app}=-\varepsilon^2u_1-2\varepsilon^3u_2.
\]
Its \(Z_t=\partial_s p\) equation then cancels identically,
and the remaining residual is \(O(\varepsilon_0^3P(T))\).
The same majorant for the error gives that order of error
for \(u,u_t\) on a smaller fixed \(z\)-disc.
Cauchy estimates and the physical derivative
\(\partial_s=-\varepsilon_0z^2\partial_z\)
improve the spatial error to \(O(\varepsilon_0^4P_1(T))\).
This establishes a genuine logarithmic-time development, not
a formal asymptotic substitution.

Take any fixed \(\gamma>3/d\) and evaluate at
\(t=\gamma\log s\). The preceding estimates give
\[
\begin{split}
 \partial_s\log a&=\mu\varepsilon+O(\varepsilon^2P(t)),&
 \partial_t\log a&=-\ell\mu\varepsilon/d+o(\varepsilon),\\
 \partial_s\log b_{\rm orb}&=N_1\varepsilon^2/d^2+
                                      o(\varepsilon^2),&
 \partial_t\log b_{\rm orb}&=-\ell N_1\varepsilon^2/d^3+
                                      o(\varepsilon^2).
\end{split}
\]
Also \(H\to1\). Therefore
\[
 \langle\nabla\log a,\nabla\log b_{\rm orb}\rangle_h
 =\frac{\mu N_1}{d^2}
       \left(1+\frac{\ell^2}{d^2}\right)s^{-3}
       +o(s^{-3}).
\]
If \(b>\mu\), then \(N_1>0\), and the mixed-circle curvature,
the negative of this quantity, is negative.
The estimates keep the normal development regular along the
whole segment; a complete soliton agrees with it by analytic
uniqueness. This contradicts positive sectional curvature and
proves \(b\le\mu\).

There is a second sectional-curvature test at the same points.
Let \(h=\log H\), \(v=\log b_{\rm orb}\), and
\(e_s=H^{-1}\partial_s\). The full connection formula is
\[
 K(e_s,e_\psi)=-H^{-2}
 \left(v_{ss}+v_s^2-h_s v_s+H^2h_t v_t\right).
\]
At fixed \(t\) the physical spatial derivative is
\(-\varepsilon^2\partial_\varepsilon\); it is taken before
evaluating at \(t=\gamma\log s\). The expansion gives
\[
 v_{ss}=-\frac{2N_1}{d^2}\varepsilon^3
              +O(\varepsilon^4P(t))+o(\varepsilon^3),\quad
 h_s=O(\varepsilon^3P(t)),\quad h_t=O(\varepsilon^2P(t)),
\]
and \(v_s=O(\varepsilon^2P(t))\),
\(v_t=O(\varepsilon^2P(t))+O((1+t)e^{-dt})\).
Here and below \(P\) denotes a fixed polynomial, whose value can
increase from line to line. Two Cauchy estimates on a smaller
complex \(z\)-disc make the twice-spatially-differentiated remainder
\(O(\varepsilon_0^5P(T))\); differentiating
\(-\varepsilon_0z^2\partial_z\) also has precisely this order.
Thus the connection terms above cannot contribute at order
\(\varepsilon^3\), and
\[
 K(e_s,e_\psi)=\frac{2N_1}{d^2}s^{-3}+o(s^{-3}).
\]
Positive sectional curvature now requires \(N_1\ge0\).
Together with the mixed-circle test this gives \(N_1=0\),
equivalently \(b=\mu\).

We next determine the second coefficient; this prevents the
first equality from leaving an arbitrary normal-radius mode.
The scalar identity on the axis gives
\[
 N=c+N_2\varepsilon^2+O(\varepsilon^3),\qquad
 N_2=\ell(\omega-g)-\frac{\mu^2}{2}.
\]
Now \(\beta_1=0\), \(Q_1'=0\), and
\[
 \alpha_1'=-\frac{\ell\mu}{d}\tanh(dt/2).
\]
Consequently the second normal subsystem is exactly
\[
 \beta_2''+D\beta_2'+LQ_2'=0,\qquad
 Q_2''=-2L\beta_2'-(\alpha_1')^2.
\]
Put \(x=dt\), \(P_0=\ell^2\mu^2/d^4\), and
\(Y(x)=\sinh x\,\partial_x\beta_2\).
The first equation gives \(Q_2'/d=-Y'\); substituting this
into the second gives the scalar equation
\[
 Y''-2\operatorname{csch}^2x\,Y=P_0\tanh^2(x/2),\qquad
 Y(x)=-\frac{N_2}{3d^2}x^2+O(x^4).
\]
The homogeneous functions \(\coth x\) and \(x\coth x-1\)
have Wronskian one. Variation of parameters gives
\[
\begin{split}
Y(x)={}&-\frac{N_2}{d^2}(x\coth x-1)\\
&+P_0(x\coth x-1)\int_0^x\coth r\,\tanh^2(r/2)\,dr\\
&-P_0\coth x\int_0^x(r\coth r-1)\tanh^2(r/2)\,dr.
\end{split}
\]
All integrals are regular at zero. This representation proves
\(Y=O((1+x)^2)\), \(Y'=O(1+x)\), and hence exponential,
polynomially weighted decay of \(\beta_2'\) and \(\beta_2''\).
In particular \(\beta_2\) has a finite limit \(E\).

That limit can be evaluated exactly. For
\(\mathcal D=\partial_x^2-2\operatorname{csch}^2x\),
the function \(z=\operatorname{csch}x\) satisfies \(\mathcal Dz=z\).
Twice integrating by parts yields
\[
\begin{split}
E=\int_0^\infty Yz\,dx
&=[Yz'-Y'z]_0^\infty+
      P_0\int_0^\infty\operatorname{csch}x\tanh^2(x/2)\,dx\\
&=-\frac{N_2}{d^2}+\frac{P_0}{2}.
\end{split}
\]
The boundary term at infinity is zero by the polynomial growth
of \(Y\), and at zero it is \(N_2/d^2\).
The integral is \(1/2\), by the substitution \(r=x/2\):
its integrand becomes \(\sinh r/\cosh^3r\).

For \(\gamma>4/d\), the same full connection formula and the
twice-spatially-differentiated second-order expansion give
\[
 K(e_s,e_\psi)=-6E s^{-4}+o(s^{-4}).
\]
Indeed \(v_{ss}=6\varepsilon^4\beta_2+O(\varepsilon^5P(t))\),
\(h_s v_s=O(\varepsilon^6P(t))\), and
\(h_t v_t=O(\varepsilon^5P(t))+O(\varepsilon^2P(t)e^{-dt})\).
This already requires \(E\le0\).

To get the reverse inequality, one more coefficient is needed
in the normal time derivative. No uncontrolled differentiation of
a pointwise error is used. In normalized variables the exact
second-radius equation is
\[
 \beta_{tt}+D\beta_t+LQ_t
 =-\beta_tQ_t-e^{-2h}(\beta_{ss}+\Omega\beta_s),
\]
where
\[
 \Omega=(\log w)'+F_0'+\alpha_s+\beta_s-h_s+f_s.
\]
Since \(b=\mu\), \(\beta_1=Q_1'=h_1=0\), and
\(f_1=-\alpha_1\), its third-order equation and the normal trace
equation are
\[
 \beta_3''+D\beta_3'+LQ_3'=2\ell\beta_2,\qquad
 Q_3''=-2L\beta_3'-2\alpha_1'\alpha_2'.
\]
In particular there is no additional nondecaying third-order
source from \(h_2\), \(Q_2\), or \(\alpha_2\) in the first equation:
each multiplies \(\beta_1\), \(\beta_1'\), or a spatial derivative
of order at least three. The latter terms affect order four or
higher. The second equation gives \(Q_3'=O((1+t)^2)\), and
variation of constants in the stable first equation gives
\[
 \beta_3'\longrightarrow\frac{2\ell E}{d}.
\]
These estimates follow jointly: the \((\beta_3',Q_3')\) homogeneous
propagator is bounded, and the inhomogeneous terms have polynomial
growth. The other two third-order coefficient equations have
the scalar stable operator \(\partial_t^2+D\partial_t\), with
polynomially bounded lower-order sources. Thus \(u_3\) and all
required time derivatives have polynomial growth as well.

The previous analytic majorant now applies with
\[
 u_{\rm app}=\varepsilon u_1+\varepsilon^2u_2+
                  \varepsilon^3u_3,\qquad
 Z_{\rm app}=-\varepsilon^2u_1-2\varepsilon^3u_2-
                  3\varepsilon^4u_3.
\]
The \(Z_t=\partial_s p\) equation again cancels identically.
All other residuals are \(O(\varepsilon_0^4P(T))\), so the
same bounded-propagator estimate gives that error for \(u,u_t\)
and one extra power for its spatial derivative, uniformly on
\(t\le\gamma\log(1/\varepsilon_0)\).
Consequently, on \(t=\gamma\log s\),
\[
 v_s=-2E\varepsilon^3+o(\varepsilon^3),\qquad
 v_t=\frac{2\ell E}{d}\varepsilon^3+o(\varepsilon^3).
\]
The previously obtained first-radius derivatives then give
\[
 K(e_\phi,e_\psi)
 =2\mu E\left(1+\frac{\ell^2}{d^2}\right)s^{-4}
       +o(s^{-4}).
\]
Positive sectional curvature requires \(E\ge0\). Hence \(E=0\),
which is exactly \(N_2=\ell^2\mu^2/(2d^2)\), and substitution
of the axis scalar identity gives the two necessary equalities
\[
 b=\mu,\qquad
 g=\omega-\frac{\mu^2}{2\ell}
                  \left(1+\frac{\ell^2}{d^2}\right).
\]

These equalities do not merely restrict the next free coefficients.
Together they force the meridional curvature to have the wrong sign.
The following fourth-order calculation makes that obstruction
independent of all remaining axis coefficients.
Write
\[
 \log W(\varepsilon)=\omega\varepsilon+\zeta\varepsilon^2
                         +O(\varepsilon^3),\qquad
 F_0'=\ell-\mu\varepsilon+g\varepsilon^2+j\varepsilon^3
                         +O(\varepsilon^4).
\]
Set
\[
 A_*=-\frac{\ell\mu}{d},\qquad
 a_0=\frac{2\ell\mu\log2}{d^2},\qquad
 \eta=-\frac{\mu^2}{d}.
\]
The already solved equations give the following polynomial parts
at large \(t\), each with an exponentially decaying, polynomially
weighted remainder (also after the derivatives used below):
\[
\begin{split}
 \alpha_1&=A_*t+a_0+O(e^{-dt}),&
 h_2&=\eta t+\frac\mu\ell a_0+O(e^{-dt}),\\
 \alpha_2&=a_{22}t^2+a_{21}t+a_{20}+O(P(t)e^{-dt}),&
 Q_2&=q_{22}t^2+q_{21}t+q_{20}+O(P(t)e^{-dt}),\\
 h_3&=\eta_3t+h_{30}+O(e^{-dt}),
\end{split}
\]
where
\[
\begin{split}
 a_{22}&=-\frac{\ell^2\mu}{2d^2},&
 a_{21}&=\frac{\mu+\ell\omega+\ell a_0+\ell^2\mu/d^2}{d},\\
 q_{22}&=-\frac{A_*^2}{2},&
 q_{21}&=\frac{\ell^2\mu^2(1+2\log2)}{d^3},\\
 \eta_3&=-\frac{2((1-\mu)\omega-g)}d.
\end{split}
\]
Here \(h_2=(\mu/\ell)\alpha_1\) holds exactly, since both
satisfy the same stable scalar equation with proportional constant
sources and regular zero data. The formula for \(q_{21}\) also
has an exact derivation, rather than a numerical fit. In the
variation-of-parameters expression for \(Y\),
\[
 \int_0^x\coth r\tanh^2(r/2)\,dr
 =-\tfrac12\tanh^2(x/2)+2\log\cosh(x/2)
 =x-\tfrac12-2\log2+O(e^{-x}).
\]
Integrating the other term by parts and using
\(N_2/d^2=P_0/2\) gives
\[
 Y=\tfrac12P_0x^2-P_0(1+2\log2)x+O(1).
\]
The relation \(Q_2'/d=-Y'\) proves the displayed \(q_{21}\).

To obtain the equation for \(h_4\), use the exact first diagonal
equation in the form
\[
 h_{tt}+Dh_t=-h_tQ_t-e^{-2h}\mathcal S,
\]
\[
\mathcal S=(\log w)''+F_0''+(Q-h)_{ss}
 +((\log w)'+\alpha_s)^2+\beta_s^2
 -h_s((\log w)'+F_0'+(Q-h)_s).
\]
Its coefficients are
\[
\begin{split}
 \mathcal S_2&=\mu^2,\\
 \mathcal S_3&=2((1-\mu)\omega-g)-2\mu\alpha_1+2\ell h_2
              =2((1-\mu)\omega-g),\\
 \mathcal S_4&=(6-4\mu)\zeta-3j+6Q_2-6h_2
      +(\omega+\alpha_1)^2-4\mu\alpha_2+3\ell h_3.
\end{split}
\]
This also verifies the constant-source equation used for \(h_3\).
Thus
\[
 h_4''+Dh_4'=-h_2'Q_2'-\mathcal S_4+2\mu^2h_2.
\]
The quadratic coefficient on its right side vanishes:
\(6q_{22}+A_*^2-4\mu a_{22}=0\).
The linear coefficient is
\[
\begin{split}
 C_t={}&\eta A_*^2-6q_{21}+6\eta
       -2(\omega+a_0)A_*+4\mu a_{21}-3\ell\eta_3+2\mu^2\eta\\
 ={}&\frac{\mu^2}{d^3}
          \big((d^2+\ell^2)-\mu^2(2d^2+\ell^2)\big).
\end{split}
\]
The last identity substitutes the forced value of \(g\).
The coefficients \(\zeta,j\) affect only the constant term.
The full source is \(C_tt+C_0+O(P(t)e^{-dt})\), so the
stable scalar equation yields \(h_4''\to C_t/d\).

The use of this coefficient at \(t=\gamma\log s\) again requires
a uniform remainder, now through four orders. The coefficient
equations are triangular at every fixed order: each physical
\(s\)-derivative raises the power of \(\varepsilon\), the two
\(h,\alpha\) homogeneous operators are
\(\partial_t^2+D\partial_t\), and the
\((\beta_t,Q_t)\) homogeneous propagator is bounded. Their
fourth-order sources contain only the first three coefficients,
their derivatives, and bounded cigar coefficients. They are
polynomially bounded. This proves polynomial bounds for \(u_4\)
and all time derivatives needed here. Use
\[
 u_{\rm app}=\sum_{r=1}^4\varepsilon^r u_r,\qquad
 Z_{\rm app}=-\sum_{r=1}^4r\varepsilon^{r+1}u_r
\]
in the same first-order system and analytic majorant.
Its residual is \(O(\varepsilon_0^5P(T))\), with exact cancellation
in \(Z_t=\partial_s p\). Hence the remainders in \(u,u_t\)
are \(O(\varepsilon_0^5P(T))\) on a fixed smaller \(z\)-disc.
Cauchy estimates control the required spatial derivatives there.
Substitution into the exact diagonal equations then gives the
same \(O(\varepsilon_0^5P(T))\) error for \(u_{tt}\);
this does not require differentiating an unquantified time error.

Choose any fixed \(\gamma>5/d\). The exponential tails of
\(h_2'',h_3''\) at \(t=\gamma\log s\) are then smaller
than order \(s^{-4}\). Since the meridional metric is
\(dt^2+H^2ds^2\), its sectional curvature is exactly
\[
\begin{split}
 K(e_s,e_t)&=-H_{tt}/H=-h_{tt}-h_t^2\\
 &=-\left(\frac{C_t}{d}+\eta^2\right)s^{-4}+o(s^{-4})\\
 &=-\frac{\mu^2(1-\mu^2)(d^2+\ell^2)}{d^4}s^{-4}+o(s^{-4}).
\end{split}
\]
For \(0<\mu<1\) the coefficient is strictly negative. The
regular normal development agrees with any complete ambient
soliton by analytic uniqueness, so this contradicts positive
sectional curvature and proves the asserted exclusion.
\(\square\)

**Consequence for the entire diagnostic family in Example 6.5.**
None of those axis data, for any of their admissible parameters,
comes from a complete positively curved steady soliton.
Indeed their analytic expansion has
\(\mu=1-2a\), \(b=B\), \(\omega=g=0\), and
\(\ell=B\pi/2\), with \(c>0\). If \(B\ne\mu\), one of the
two order-\(s^{-3}\) sectional curvatures above is negative.
If \(B=\mu\), then \(N_2=-\mu^2/2\), so
\[
 E=\frac{\mu^2}{2d^2}\left(1+\frac{\ell^2}{d^2}\right)>0,
 \qquad K(e_s,e_\psi)=-6Es^{-4}+o(s^{-4})<0.
\]
This last case uses only the controlled second-order development;
the third-order calculation is needed for the general exact
coefficient theorem, not for excluding this family.

The original examples still show that axis positivity and scalar
identities alone are insufficient. The stronger ambient curvature
test now excludes every analytic sublinear power-law axis with
\(c>0\), including all these examples. Actual LP axes have not
been shown to have the expansion assumed here; a general
sublinear axis or the branch \(c=0\) is not excluded by this
argument. The strict boundary-orbit gap remains open.

### Corollary 6.7a (finite holomorphic asymptotics suffice)

The exclusion in Proposition 6.7 does not require a convergent power series
at infinity. It still holds under the following finite, uniform hypothesis.
For some fixed `0<rho<1` and `0<eta<=1`, the positive real axis functions
`log w` and `f=F_0'` extend holomorphically to every disc
`|z-s|<rho*s`, for all sufficiently large real `s`, and there have expansions

\[
 \log w(z)=\log w_0+\mu\log z+
                  \sum_{r=1}^4 w_rz^{-r}+O(s^{-4-\eta}),
 \qquad
 f(z)=\ell-bz^{-1}+\sum_{r=2}^4 f_rz^{-r}
                  +O(s^{-4-\eta}).
\]

The constants in both errors are independent of `s` and `z` in those
discs; the coefficients are fixed real numbers. Retain `w_0>0`,
`0<mu<1`, `ell>0`, and positive limiting normal sectional curvature.
Then no complete positive-sectional-curvature extension exists.

**Proof.** On each smaller concentric disc, Cauchy's inequalities give
`O(s^{-4-eta-k})` for the `k`-th spatial derivative of either remainder.
Thus the fourth-order formal coefficients in Proposition 6.7 are unchanged,
with `omega=w_1`, `zeta=w_2`, `g=f_2`, and `j=f_3`. All normal-axis data
are obtained from these functions and their spatial derivatives by the
axis equations; their normalized errors are `O(s^{-4-eta})` in the same
holomorphic norm. No derivatives of an error on the real line alone are
being inferred here.

For a large fixed `s_0`, write `epsilon_0=1/s_0` and use
`xi=s_0/z` as the spatial variable. A fixed disc about `xi=1` lies in
the assumed holomorphic domain, independently of `s_0`; the physical
derivative is `-epsilon_0 xi^2 partial_xi`. The parameter-dependent
regular-singular construction in the proof of Proposition 6.7 therefore
has the same uniform initial interval. Its coefficientwise analytic
majorant also applies to differences of the initial axis data: subtract
the two recurrences, use the same inverse bounds, and majorize each
analytic difference by its first derivative times the data difference
on a slightly smaller fixed polydisc. This gives an
`O(epsilon_0^{4+eta})` error in the normalized initial development.

The fourth-order approximate solution has equation residual
`O(epsilon_0^5 P(T))`. In the logarithmic-time difference estimate from
Proposition 6.7, retain also the preceding initial error. Its bounded
linear propagator and analytic radius estimate give

\[
 \|u-u_{\rm app}\|+\|u_t-(u_{\rm app})_t\|
 \le C\bigl(\epsilon_0^{4+\eta}
                +\epsilon_0^5\bigr)P_1(T),
 \qquad 0\le t\le T=\gamma\log s_0.
\]

The same fixed-disc Cauchy estimates give an extra factor `epsilon_0`
per spatial derivative, and the exact diagonal equations bound the
second normal derivative. The right-hand side is `o(s_0^{-4})` for every
fixed `gamma`, since `eta>0`. The estimates at orders two and three
likewise have the errors needed for the earlier two sign comparisons.
They force `b=mu` and
`g=omega-mu^2(1+ell^2/d^2)/(2ell)`, exactly as before.
The fourth-order meridional curvature is consequently

\[
 K(e_s,e_t)=-\frac{\mu^2(1-\mu^2)(d^2+\ell^2)}{d^4}s^{-4}
                         +o(s^{-4})<0.
\]

Analytic uniqueness identifies this regular development with that of an
actual complete soliton, giving the same contradiction. \(\square\)

This is a finite-asymptotic version of the same obstruction, not a proof
of the needed expansions for the LP examples. For instance, exponentially
small holomorphic perturbations satisfy its remainder condition even
when they have an essential singularity in `1/s` at zero. Real-axis
asymptotics without the uniform complex-disc bounds are insufficient for
this argument; the general sublinear and zero-limit branches remain open.
The exact algebra, including the connection square in `-H_tt/H`, is
reconstructed by `scripts/check_sublinear_fourth_order.py`. That script
does not replace the holomorphic estimates.

### Proposition 6.8 (global transverse Jacobi budget)

Let the complete positively curved steady soliton carry the standard
\(O(2)\times O(2)\) action, meaning an action equivariantly
diffeomorphic to the linear action on
\(\mathbb R^2\oplus\mathbb R^2\). Let \(F\) be its invariant
potential. From \(p\in\Sigma_p\setminus\{o\}\), launch any unit
geodesic \(\gamma\) normal to \(\Sigma_p\), and set
\(w=\mathcal A(p)>0\). Then \(\gamma\) first meets \(\Sigma_q\)
at a finite time \(T>w\), away from \(o\). Writing
\[
 a(t)=\mathcal A(\gamma(t)),\qquad
 \kappa(t)=\sec(\gamma'(t),K_\phi(t)/a(t)),\quad 0\le t<T,
\]
where \(K_\phi\) is the period-\(2\pi\) first-circle generator,
there is a number \(0<\nu<1\) such that
\[
 \int_0^T\kappa(t)a(t)\,dt=\nu,\qquad
 \int_0^L\kappa(t)\,dt<\frac1{w-L}\quad(0<L<w).
\]
In particular,
\[
 \int_0^{w/2}\kappa(t)\,dt<\frac2w.
\]
No global gradient bound for \(\mathcal A\), injectivity of the
normal exponential map, or curvature-scale compactness is needed.

**Proof.** Choose a reflection in each factor fixing the initial
point and the selected initial tangent vector. Their common fixed
surface is a smooth totally geodesic meridional plane.
Geodesic uniqueness keeps \(\gamma\) in that plane.
On the principal part, the reflections eliminate angular cross
terms and the metric is
\(h+\mathcal A^2d\phi^2+\mathcal B^2d\psi^2\).
Along a meridional geodesic the unit angular field
\(K_\phi/\mathcal A\) is parallel. The Killing Jacobi equation
therefore gives the scalar equation
\[
 a''=-\kappa a.
\]
The reflection fixing \(\Sigma_p\) pointwise reverses the initial
normal direction and preserves \(\mathcal A\); thus
\(a(0)=w\) and \(a'(0)=0\).

The standard linear action has no exceptional orbit strata:
the nonprincipal set consists exactly of the two coordinate
fixed planes. Before encountering \(\Sigma_q\), the field
\(K_\phi\) is nonzero. Both \(a\) and the scalar Jacobi equation
extend smoothly across any intermediate crossings of \(\Sigma_p\).
Such a crossing cannot be tangent, by total geodesicity and
uniqueness. Focal points of the initial normal family do not affect
this equation for the individual geodesic and its Killing field.

Since \(\kappa>0\), the derivative \(a'\) strictly decreases
as long as \(a>0\). If no zero occurred, fixing any \(t_0>0\)
would give
\[
 a(t)\le a(t_0)+(t-t_0)a'(t_0),\qquad t\ge t_0,
\]
with \(a'(t_0)<0\), a contradiction for sufficiently large \(t\).
Completeness prevents finite-time escape. Thus there is a finite
first zero \(T\), whose endpoint lies in \(\Sigma_q\).

Let \(\nu\) be the magnitude of the component of \(\gamma'(T)\)
normal to \(\Sigma_q\). Smoothness of the standard polar action,
with unit rotation weight and period \(2\pi\), gives
\(a'(T-)=-\nu\), \(0\le\nu\le1\).
If \(\nu=0\), geodesic uniqueness at this totally geodesic
fixed plane would force \(\gamma\) to have lain in it, which is
impossible. Hence \(\nu>0\).

Invariance makes \(\nabla F\) tangent to both fixed planes.
Consequently
\[
 (F\circ\gamma)'(0)=0,\qquad
 (F\circ\gamma)'(T)
   =\int_0^T\operatorname{Ric}(\gamma',\gamma')\,dt>0.
\]
If \(\nu=1\), the terminal tangent would be normal to
\(\Sigma_q\), forcing the latter derivative to vanish.
This proves \(\nu<1\). The same argument excludes a hit at
the common vertex, where the invariant potential has zero
gradient. An invariant potential can also be obtained by
averaging over the compact group: its Hessian remains Ricci,
and positive Ricci curvature rules out a nonzero parallel
difference of potential gradients.

Integration of the Jacobi equation, with the smooth terminal
slope just established, gives
\(\int_0^T\kappa a=\nu\).
Strict decrease of \(a'\) from zero to \(-\nu\) also gives
\[
 a(t)>w-\nu t\ge w-t\quad(0<t<T),\qquad
 w=\int_0^T(-a')\,dt<\nu T<T.
\]
Thus \(T>w\). For \(0<L<w\), integration on \([0,L]\) yields
\[
 (w-L)\int_0^L\kappa\,dt
 \le\int_0^L\kappa a\,dt<\nu<1,
\]
as required. \(\square\)

This first collision need not minimize distance between the
fixed planes; in fact its tangent is not normal to the second
plane. The conclusion does not provide global normal coordinates
or solve the second-axis soliton matching problem.

**Conditional consequence.** Suppose an actual axis satisfies
\(w(s)\sim w_0s^\mu\). If fixed constants \(\delta,k_0>0\) give
the additional transverse estimate
\[
 \kappa_s(t)\ge k_0/s
 \quad\text{for }0\le t\le\delta\sqrt{s}
\]
for all sufficiently large \(s\), then \(\mu\le1/2\).
Indeed, if \(\mu>1/2\), that whole interval is contained in
\([0,w(s)/2]\), and
\[
 \frac{k_0\delta}{\sqrt{s}}
 \le\int_0^{\delta\sqrt{s}}\kappa_s(t)\,dt
 <\frac2{w(s)},
\]
which is impossible as \(s\to\infty\).
The displayed transverse lower bound has **not** been proved
for the LP family. An asymptotic value of \(\kappa_s(0)\)
on the axis alone does not imply it. The proposition supplies
a global integral constraint on an expanding transverse interval,
not a resolution of the zero-limit case or the strict gap.

## 7. The remaining analytic obstacle

For \(p=q=2\), the original boundary question is now the precise strict
curvature-gap problem
\[
 \int_{\Sigma_p}K_{\Sigma_p}\,dA<2\pi.
\]
The estimates above prove only the non-strict upper bound \(2\pi\).
Strict positivity of Gaussian curvature does not by itself make that
bound strict: the cigar has positive curvature and total curvature
\(2\pi\). The product examples in Section 5 also preserve any prescribed
tip ratio while attaining this value. They do not have strictly positive
ambient curvature operator, so exploiting that extra hypothesis is the
remaining task.

### 7.1. The exact axis equation

On \(\Sigma_p\) in dimension four, let \(e_s,e_\phi\) be its unit radial
and angular vectors, and let \(\nu\) be any unit normal vector. Write
\[
 k=K(e_s,e_\phi)=-w''/w,\qquad
 U=K(e_s,\nu),\quad V=K(e_\phi,\nu),\quad
 N=K(\nu_1,\nu_2).
\]
The \(O(2)\) normal isotropy makes \(U,V\) independent of the choice of
\(\nu\). Restricting the tangential steady equations gives
\[
 F''=k+2U,\qquad F'\frac{w'}w=k+2V,
 \quad\hbox{hence}\quad w''+F'w'=2wV.
 \tag{20}
\]
With \(F(0)=0\), integration yields
\[
 w'(s)=e^{-F(s)}\left(1+2\int_0^s e^{F(t)}w(t)V(t)\,dt\right).
 \tag{21}
\]
Strict positivity of \(V\), even divergence of the integral, does not
make the limit of this quotient positive.

### Proposition 7.1 (positive-curvature axis jets do not force a gap)

For every \(0<\theta<1\), there are smooth standard
\(O(2)\times O(2)\)-invariant metric and potential jets along a complete
fixed surface \(\Sigma_p\) such that:

1. the curvature operator is strictly positive along \(\Sigma_p\);
2. the full tensor equation \(\operatorname{Ric}=\nabla^2F\) holds at every
   point of \(\Sigma_p\);
3. the scalar identity \(R+|\nabla F|^2=\mathrm{constant}\) holds there;
4. \(\lambda_p/\lambda_q=\theta\) at the vertex;
5. nevertheless \(c_p=0\) and \(\int_{\Sigma_p}K\,dA=2\pi\).

These jets are realized by a smooth positively curved metric on an open
neighborhood of \(\Sigma_p\). They are **not** asserted to solve the
steady equation off \(\Sigma_p\), nor to extend to a complete ambient
soliton.

**Proof.** Put
\[
 w(s)=\tanh s,\quad F_0(s)=3\log\cosh s,\quad S(s)=\operatorname{sech}^2s,
\]
\[
 k=2S,\quad U=V=S/2,\quad
 N=2/\theta-1/2+S/2.
 \tag{22}
\]
Every displayed curvature is strictly positive. Direct differentiation
gives both equations (20). Also
\[
 R:=2k+4U+4V+2N=9S+4/\theta-1,
 \quad R+(F_0')^2=8+4/\theta,
\]
and \(k(0)/N(0)=\theta\). Since \(w'=S\to0\), its total Gaussian
curvature is \(2\pi\). Formula (21) here reads
\[
 1+2\int_0^s e^{F_0}wV\,dt=\cosh s,
 \qquad e^{F_0(s)}=\cosh^3s,
\]
so even the diverging numerator fails to give a positive limit.

To realize these data geometrically, use polar coordinates \((t,\psi)\)
normal to \(\Sigma_p\), and set
\[
 H(s,t)=1-\tfrac12U(s)t^2,\qquad
 B(s,t)=t-\tfrac16N(s)t^3,
\]
\[
 \widehat g=dt^2+H(s,t)^2\bigl(ds^2+w(s)^2d\phi^2\bigr)
                       +B(s,t)^2d\psi^2,
\]
\[
 \widehat F(s,t)=F_0(s)+\tfrac12(2U(s)+N(s))t^2.
 \tag{23}
\]
These tensors are smooth at \(t=0\): \(H\) is even, \(B\) is odd with
\(B_t(s,0)=1\). They are also smooth at \(s=0\), including the common
vertex: \(U,N,F_0\) are smooth even radial functions, while \(w\) is odd
with \(w'(0)=1\). Choose an open neighborhood where \(H>0\) and \(B>0\)
for \(t>0\).

At \(t=0\), the curvature operator in the natural two-plane basis is
diagonal with eigenvalues
\[
 k,U,U,U,U,N.
\]
Indeed, \(-H_{tt}/H=U\), and the first circle radius \(a=wH\) also gives
\(-a_{tt}/a=U\), at \(t=0\). Moreover \(-B_{tt}/B\to N\), and the tangent
two-plane has curvature \(k\). Normal rotation symmetry supplies the
paired repetitions. The tangential Ricci components are \(k+2U\), matching the
tangential Hessian by (20) and \(U=V\). Each normal Ricci component is
\(2U+N\), matching \(\widehat F_{tt}(s,0)\); cross terms vanish by the
reflections. Thus the full steady tensor equation holds on the axis.
The scalar identity follows from the displayed calculation. Positivity
of the curvature operator persists on some open neighborhood of the
entire axis by continuity, with neighborhood width allowed to vary.
\(\square\)

### Proposition 7.2 (formal compatibility to every transverse order)

For each \(0<\theta<1\), the axis data in Proposition 7.1 extend
uniquely in the normal-distance gauge to formal series
\[
\begin{aligned}
 H&=1+\sum_{j\ge1}h_{2j}(s)t^{2j},&
 a&=w(s)\left(1+\sum_{j\ge1}a_{2j}(s)t^{2j}\right),\\
 b&=t\left(1+\sum_{j\ge1}b_{2j}(s)t^{2j}\right),&
 F&=F_0(s)+\sum_{j\ge1}f_{2j}(s)t^{2j},
\end{aligned}
\]
such that
\[
 g=dt^2+H^2ds^2+a^2d\phi^2+b^2d\psi^2,
 \qquad \operatorname{Ric}_g-\nabla^2F=0
\]
as a formal series. The initial coefficients are
\[
 h_2=a_2=-U/2,\qquad b_2=-N/6,\qquad f_2=(2U+N)/2.
\]
All coefficients are smooth along the complete fixed surface, including
its vertex. For every finite prescribed order, a sufficiently long
truncation gives a smooth metric with positive curvature operator on
some neighborhood of that surface and with soliton residual vanishing
to that order in the normal distance.

This is a statement about formal series and finite truncations.
**Convergence, an exact soliton on an open set, and a complete ambient
soliton are not asserted.**

**Proof.** Write \(E=\operatorname{Ric}-\nabla^2F\).
First work at \(s>0\). Reflections in \(\phi,\psi\) leave only the four
diagonal components and \(E_{st}\) potentially nonzero. The diagonal
components \(E_{ss}/H^2,E_{\phi\phi}/a^2,E_{\psi\psi}/b^2,E_{tt}\)
are even formal series, whereas \(E_{st}\) is odd.

Suppose all coefficients below an even order \(n\ge4\) have been
constructed. Add
\[
 \delta H=h_nt^n,\quad \delta a=w a_nt^n,\quad
 \delta b=b_nt^{n+1},\quad \delta F=f_nt^n.
\]
At order \(t^{n-2}\), the contributions of these new coefficients to the
four normalized diagonal components are
\[
\begin{array}{c|l}
 ss/H^2 &-n^2h_n\\
 \phi\phi/a^2&-n^2a_n\\
 tt&-n(n-1)(h_n+a_n+f_n)-n(n+1)b_n\\
 \psi\psi/b^2&-n(h_n+a_n+f_n)-n(n+1)b_n .
\end{array}
\tag{24}
\]
All other terms at that order depend only on earlier coefficients and
their \(s\)-derivatives. To see (24), the second \(t\)-derivatives
contribute \(n(n-1)\) for \(H,a,F\) and \(n(n+1)\) for \(b/t\).
The leading \(b_t/b=1/t+O(t)\) contributes the additional \(n\) in
the first two rows and the terms in the last row. Derivatives solely
in \(s\) cannot lower the new power of \(t\).

The first two rows determine \(h_n,a_n\). The last two then determine
\(b_n,f_n\): their coefficient matrix has determinant of absolute
value \(n^2(n+1)(n-2)>0\). Thus all diagonal equations have a unique
formal solution extending the specified order-two data.

It remains to prove the mixed equation, not merely the diagonal ones.
Set \(C=R+|\nabla F|^2\). The contracted Bianchi identity gives
\[
 \operatorname{div}E-d(\operatorname{tr}E)+E(\nabla F,\cdot)
       =-\tfrac12dC.
 \tag{25}
\]
For the formal diagonal solution, only \(e=E_{st}\) remains, so
\[
 C_s=-2(e_t+A e),\qquad
 C_t=-2H^{-2}(e_s+B e),
 \tag{26}
\]
where
\[
 A=\frac{H_t}{H}+\frac{a_t}{a}+\frac{b_t}{b}+F_t
       =t^{-1}+O(t),\qquad
 B=-\frac{H_s}{H}+\frac{a_s}{a}+\frac{b_s}{b}+F_s=O(1).
\]
Here both order estimates are as \(t\to0\) for fixed \(s>0\).
These identities follow by taking the divergence of the tensor whose
only nonzero entries are \(E_{st}=E_{ts}=e\).
The value of \(C\) on \(t=0\) is the constant \(8+4/\theta\) from
Proposition 7.1; higher coefficients do not change it.
Consequently (26) forces the coefficient of \(t\) in \(e\) to vanish.

If \(e_k(s)t^k\), with odd \(k\ge3\), were its first nonzero term,
commuting the two derivatives of \(C\) in (26) would give
\[
 \partial_t(e_t+A e)
       =\partial_s\!\left[H^{-2}(e_s+B e)\right].
\]
The left side has coefficient \((k+1)(k-1)e_k\) at order \(t^{k-2}\);
the right side starts at order \(t^k\). Thus \(e_k=0\), a
contradiction. All mixed coefficients vanish. This proves the full
formal equation, not just its diagonal restriction.

For smoothness at \(s=0\), induct on the truncation order. Initially the
coefficients are even smooth functions of \(s\), with \(h_2(0)=a_2(0)\).
A smooth rotational metric in the first plane has \(H\) and \(a/s\)
even and \(H(0,t)=(a/s)(0,t)\): their squared difference is divisible
by \(s^2\), which is exactly the Cartesian smoothness condition.
Assume it holds for the previous truncation. Its residual \(E\) is a
smooth invariant tensor. Hence its normalized radial and angular
tangential entries are smooth even functions and agree at \(s=0\).
Their next coefficient sources therefore agree there. The identical
first two coefficients \(-n^2\) in (24) enforce
\(h_n(0)=a_n(0)\), while preserving even smoothness. The normal diagonal
sources are also smooth even functions of \(s\), so the remaining
linear system preserves that property for \(b_n,f_n\).
This establishes the induction. Evenness in \(t\), together with
\(b/t=1+O(t^2)\), gives smoothness in the second Cartesian plane as well.

Finally, any finite truncation is smooth and retains exactly the
axis curvature operator of Proposition 7.1. Its strict positivity
persists on an open neighborhood, whose width may depend on \(s\)
and on the truncation order. Since the formal residual is zero,
taking enough coefficients makes each finite truncation satisfy the
equation to any specified finite transverse order. \(\square\)

### Lemma 7.3 (analytic regular-singular Cauchy problem)

Let \(A\) be a constant finite-dimensional matrix whose eigenvalues have
positive real parts. If \(G(t,s,V,Z)\) is analytic near
\((0,s_0,0,0)\), then
\[
 (t\partial_t+A)V=tG(t,s,V,\partial_sV),\qquad V(0,s)=0
 \tag{27}
\]
has a unique local jointly analytic solution.

**Proof.** The coefficient equation for \(V=\sum_{n\ge1}V_n(s)t^n\) is
\[
 (nI+A)V_n=[t^{n-1}]G(t,s,V,\partial_sV).
\]
For some constant \(C\), the induced maximum norm satisfies
\(\|(nI+A)^{-1}\|\le C/n\) for all integers \(n\ge1\).
This follows by a Neumann expansion for large \(n\), then increasing
\(C\) for the finitely many remaining invertible matrices.

Expand \(G\) in powers of \(t,\xi=s-s_0,V,Z\), replace its coefficients
by their absolute values, sum over output components, and replace all
components of \(V,Z\) by scalar variables \(v,z\).
This gives an analytic scalar majorant \(\widehat G(t,\xi,v,z)\).
The ordinary analytic Cauchy--Kowalevski theorem solves
\[
 \partial_tv=C\widehat G(t,\xi,v,\partial_\xi v),\qquad v(0,\xi)=0.
\]
Its formal coefficients are nonnegative. Induction in \(n\) shows that
\(v_n(\xi)\) coefficientwise majorizes every component of \(V_n\):
differentiation in \(\xi\) preserves majorization, and the preceding
matrix bound gives the required factor \(C/n\).
Thus the vector formal series converges in a smaller bidisc.
Its recurrence proves uniqueness. \(\square\)

### Proposition 7.4 (convergent local solitons with the collapsed axis)

For every \(0<\theta<1\), the formal series in Proposition 7.2 converge
to an analytic metric and potential on an invariant open neighborhood
of the entire complete fixed surface. They satisfy
\(\operatorname{Ric}=\nabla^2F\) exactly and, after shrinking that
neighborhood, have strictly positive curvature operator.
The fixed surface still has \(w(s)=\tanh s\), total curvature \(2\pi\),
and tip ratio \(\lambda_p/\lambda_q=\theta\).

The ambient neighborhood is **not asserted complete**.
Proposition 7.5 below obtains a uniform neighborhood width, but a
complete positive-curvature extension remains open.

**Proof.** First fix \(s_0>0\). Put \(D=t\partial_t\) and
\[
 u=(x,y,z,f)=\left(\log H,\log(a/w),\log(b/t),F-F_0\right),
 \quad p=Du,\quad W=t\partial_su,\quad Y=(u,p,W).
\]
Treat \(W\) initially as independent and adjoin \(Du=p\) and
\(DW=W+t\partial_sp\). Multiplying the four normalized diagonal
equations by \(t^2\) gives an invertible constant coefficient matrix
for their second \(D\)-derivatives: the \(ss,\phi\phi,\psi\psi\)
equations contain respectively \(-Dp_x,-Dp_y,-Dp_z\), while the
\(tt\) equation contains \(-D(p_x+p_y+p_z+p_f)\).
Their regularized spatial derivatives are
\[
 tH_s/H=W_x,\quad ta_s/a=tw'/w+W_y,\quad
 tb_s/b=W_z,\quad tF_s=tF_0'+W_f.
\]
Second spatial derivatives give \(t\partial_sW\), quadratic products,
and analytic background terms with a factor \(t^2\).
Since \(w(s_0)\ne0\), solving for \(Dp\) gives the analytic system
\[
 DY=MY+Q(s,Y)+tG_0(t,s,Y,\partial_sY),
 \qquad Q(s,Y)=O(|Y|^2).
 \tag{28}
\]
The constant indicial matrix is explicitly
\[
 M(u,p,W)=(p,Bp,W),\qquad
 Bp=(0,0,-p_x-p_y-p_z-p_f,\ 2p_x+2p_y+2p_f).
\]
For example, the time-linearized diagonal equations are
\(D^2x=D^2y=0\),
\(D^2z+Dz+Dx+Dy+Df=0\), and
\((D^2-D)(x+y+f)+(D^2+D)z=0\).
Thus the spectrum of \(M\) is contained in \(\{-1,0,1,2\}\).

Let \(\bar u=u_2(s)t^2+u_4(s)t^4\) be the compatible logarithmic jet
from Proposition 7.2 and set \(P=(\bar u,D\bar u,t\partial_s\bar u)\).
Its last block includes a term of degree five.
Its residual in (28) is \(O(t^6)\), and its first and last block
equations hold exactly. Substituting \(Y=P+t^4V\) therefore yields
\[
 (D+4I-M)V=tG_1(t,s,V,\partial_sV),
\]
with \(G_1\) analytic. Indeed, \(P=O(t^2)\) makes the linearized
quadratic term \(O(t^2)V\) after division by \(t^4\); its nonlinear
remainder and the source residual also have positive powers of \(t\).
Spatial derivatives retain their factor \(t\).
The eigenvalues of \(4I-M\) lie in \(\{2,3,4,5\}\), so Lemma 7.3 applies.

The derivative constraints are preserved: \(p=Du\), and
\(Z=W-t\partial_su\) satisfies \(DZ=Z\) with \(Z=O(t^4)\), hence \(Z=0\).
Thus the solution satisfies the actual diagonal equations. Uniqueness
identifies its Taylor series with Proposition 7.2's series. It is even
in \(t\), so the metric extends analytically across the normal polar
axis. The mixed residual has every Taylor coefficient zero by
Proposition 7.2 and is analytic; it therefore vanishes identically.
This proves convergence and the full soliton equation near \(s_0>0\).

To include the common vertex, choose \(0<r<R<\pi/2\).
The background functions \(\tanh s,\operatorname{sech}^2s\), and a
branch of \(\log\cosh s\) are holomorphic on \(|s|<R\).
The formal recurrence uses differentiation, algebra, division by
nonzero constants, and formal inverses with leading terms \(1\) or
\(w(s)\). Every coefficient is consequently meromorphic on this same
disk, with only a possible finite-order pole at \(s=0\).
Proposition 7.2's real-axis smoothness rules out that pole.
The logarithmic coefficients are therefore holomorphic on a common
disk independent of their order.

The preceding convergence proof also applies at each complex point
on \(|s|=r\). Shrinking these local bidiscs to relatively compact products, a finite
covering of the circle supplies a common positive \(t\)-radius \(\varepsilon\) and finite bound \(C_0\).
Cauchy's coefficient estimates give
\[
 \sup_{|s|=r}|u_n(s)|\le C_0\varepsilon^{-n}.
\]
Maximum modulus extends this bound to \(|s|\le r\); the series
therefore converges normally on smaller bidiscs across \(s=0\).
The formal Cartesian regularity identities pass to the sum:
the functions are even in \(s,t\),
\(H(0,t)=(a/s)(0,t)\), and \((b/t)(s,0)=1\).
In particular \(H^2-(a/s)^2\) is divisible by \(s^2\), and
\(1-(b/t)^2\) by \(t^2\), as analytic functions.
These identities give an analytic metric in Cartesian coordinates
on both rotational planes. The soliton equation extends across
their axes by continuity.

Finally the unique analytic germs agree on sufficiently small
overlaps, since they have the same convergent normal Taylor series.
They glue to a neighborhood of all \(s\ge0\).
Their axis curvature operator is exactly that of Proposition 7.1,
so positivity holds after shrinking this neighborhood pointwise.
Its induced fixed-surface metric is unchanged and complete.
No step establishes ambient completeness. \(\square\)

### Proposition 7.5 (uniform positive-curvature tube and its end)

For each fixed \(0<\theta<1\), the solution in Proposition 7.4 is
defined with positive curvature operator on a tube
\[
 s\ge0,\qquad 0\le t<\varepsilon
\]
for some constant \(\varepsilon>0\) depending on \(\theta\), with the
usual smooth circle collapses at \(s=0\) and \(t=0\).
Put \(c=2/\theta-1/2\) and \(A=\sqrt{2/c}\). On any smaller fixed tube,
as \(s\to\infty\), its metric coefficients converge exponentially to
\[
 g_\infty=ds^2+d\phi^2+dt^2+B(t)^2d\psi^2,\qquad
 B(t)=A\tanh(t/A),
 \tag{29}
\]
and
\[
 F(s,t)-3s\longrightarrow
       -3\log2+2\log\cosh(t/A).
 \tag{30}
\]
The error is \(O(e^{-2s})\), also after any fixed number of
\(s\)-derivatives and Cartesian normal derivatives.
Thus the end of this tube approaches a flat cylinder times a cigar.
This proposition does not supply a complete ambient manifold;
Proposition 7.8 below gives additional finite-width end continuation.

**Proof.** Write \(\xi=e^{-2s}\). This is a parameter for the metric
coefficients in the \(ds,d\phi\) and normal frames; the metric is not
claimed to extend nondegenerately across \(\xi=0\).
The axis data satisfy
\[
 w=\frac{1-\xi}{1+\xi},\qquad
 S=\frac{4\xi}{(1+\xi)^2},\qquad
 F_0-3s=-3\log2+3\log(1+\xi).
\]
They are analytic at \(\xi=0\), and \(w(0)=1\).
In the proof of Proposition 7.4, replace each spatial derivative by
\(\partial_s=-2\xi\partial_\xi\).
The potential enters the equations through derivatives, so the
affine part \(3s\) introduces only the constant \(3\).
All coefficients and the compatible finite jet are analytic in
\(\xi\). The same indicial matrix and Lemma 7.3 therefore give
jointly analytic metric coefficients and \(F-3s\) near
\((\xi,t)=(0,0)\), with a uniform positive normal radius.
For \(\xi>0\), uniqueness identifies this solution with Proposition
7.4's solution.

At \(\xi=0\), the initial jets are \(H=1+O(t^4)\),
\(a=1+O(t^4)\),
\(b/t=1-ct^2/6+O(t^4)\), and
\(F-3s=-3\log2+ct^2/2+O(t^4)\).
The product (29), with potential (30) plus \(3s\), satisfies the exact
steady equation and these data. Indeed its cigar curvature is
\(c\,\operatorname{sech}^2(t/A)\), and the cigar potential has that
curvature times the cigar metric as Hessian.
Formal uniqueness, followed by analyticity, identifies the entire
limiting solution. Analyticity in \(\xi\) now proves the stated
exponential convergence on smaller fixed tubes. The same conclusion
holds after the stated derivatives, since
\(\partial_s=-2\xi\partial_\xi\).

Uniform positivity needs more than convergence to this product, whose
curvature operator has a five-dimensional kernel.
Use an analytic orthonormal frame along a fixed normal radial ray,
extended smoothly across \(t=0\); rotation symmetry handles the other
normal directions. Decompose bivectors into that limiting kernel and
the cigar two-plane. The curvature operator has the block form
\[
 \mathcal R(\xi,t)=
 \begin{pmatrix}
   \xi P(t)+O(\xi^2)&O(\xi)\\
   O(\xi)&c\,\operatorname{sech}^2(t/A)+O(\xi)
 \end{pmatrix}.
\]
All estimates are uniform on a sufficiently small closed \(t\)-interval.
The axis curvatures \(k=2S\), \(U=V=S/2\) show
\[
 P(0)=\operatorname{diag}(8,2,2,2,2).
\]
Thus \(P(t)\ge mI\) for some \(m>0\) after shrinking the normal interval,
and the last scalar block is bounded away from zero.
Its Schur complement is at least
\(\xi mI-C\xi^2I>0\) for sufficiently small \(\xi>0\).
This proves strict positive curvature operator uniformly in normal
radius at the end. It does not give a uniform positive curvature
lower bound as \(s\to\infty\).

The remaining finite segment of the axis is compact. Proposition
7.4 and continuity supply a common positive normal radius there,
including the vertex. Taking the minimum of that radius and the
one obtained at infinity completes the proof. \(\square\)

### Lemma 7.6 (a growing transverse mode without a sign obstruction)

Write the first coefficient of the convergent end expansion as
\(\log a=-\xi y(t)+O(\xi^2)\).
On the tube of Proposition 7.5 it satisfies
\[
 y''+D(t)y'-2y=0,\qquad
 D(t)=\frac{2}{A}\coth(2t/A),\qquad y(0)=2,\quad y'(0)=0.
 \tag{31}
\]
The regular solution of this linear ODE extends to every \(t\ge0\).
It satisfies \(y>0\), \((\log y)''>0\), and
\[
 \lim_{t\to\infty}\frac{\log y(t)}t
   =\kappa:=\frac{\sqrt{2c+8}-\sqrt{2c}}2>0.
 \tag{32}
\]
The corresponding first-order curvature block on
\(\operatorname{span}\{e_s\wedge e_\phi,e_t\wedge e_\phi\}\) is
\[
 \xi
 \begin{pmatrix}4y&-2y'\\-2y'&y''\end{pmatrix},
 \tag{33}
\]
whose coefficient matrix is positive definite at every finite \(t\).

**Proof.** Linearize the normalized \(\phi\phi\) steady equation about
the product in (29). Using \(a=1-\xi y+O(\xi^2)\) and
\(\partial_s\xi=-2\xi\) gives
\(y''+(B'/B+\overline F')y'-2y=0\), where
\(\overline F=2\log\cosh(t/A)\).
Its drift is \(D(t)\) in (31). Axis data and regularity give
\(y=2+t^2+O(t^4)\).
The regular solution exists locally by its even power series and
extends on \(t>0\) by linear ODE theory.

Set \(q=y'/y\) while \(y>0\). Then
\[
 q'=2-Dq-q^2.
\]
The positive root \(r(t)\) of \(2-Dr-r^2=0\) increases strictly,
because \(D'<0\). Near zero, \(q=t+O(t^3)\) and
\(r=2t+O(t^3)\). A first crossing is impossible: at \(q=r\) one
has \(q'=0<r'\). Nor can \(q\) cross zero, since its derivative
there is \(2\). Hence \(0<q<r\), so \(q'>0\), and \(y\) stays
positive. The increasing bounded function \(q\) tends to \(\kappa\):
a smaller limit would make \(q'\) tend to a strictly positive
constant. Integrating \(q=(\log y)'\) proves (32).

The warped curvature formula \(-\nabla_h^2a/a\) gives (33).
Its determinant is \(4y^2q'>0\) and its first entry is positive;
at zero the matrix is \(\operatorname{diag}(8,2)\).
\(\square\)

In this lemma, (31)--(33) describe the continued **linearized mode**.
Proposition 7.8 below supplies nonlinear end continuation at each
fixed finite width. This particular block
does not give a first-order curvature-loss obstruction. The mode is
unbounded, so a perturbation estimate uniform over the whole
transverse half-line cannot be inferred from this expansion.
Neither fact decides nonlinear continuation or ambient completeness.

### Lemma 7.7 (continuation along a characteristic parameter boundary)

Consider a finite-dimensional first-order analytic system
\[
 \partial_tY=G(t,\xi,Y)+\xi\,B(t,\xi,Y)\partial_\xi Y.
 \tag{34}
\]
Suppose \(Y_0(t)\) is an analytic solution of
\(Y_0'=G(t,0,Y_0)\) on a neighborhood of a compact real interval
\([t_0,T]\), and \(G,B\) are analytic on a neighborhood of its graph.
Any analytic solution germ at \((t_0,0)\) with value \(Y_0(t_0)\)
extends uniquely to a neighborhood of the whole segment
\(\{(t,0):t_0\le t\le T\}\).
The neighborhood may narrow in \(\xi\) as \(t\) increases.

**Proof.** At a continuation point \(t_*\), rescale
\(\xi=\rho z\). Equation (34) becomes
\[
 Y_t=G(t,\rho z,Y)+zB(t,\rho z,Y)Y_z.
\]
No inverse power of \(\rho\) occurs.
By shrinking \(\rho>0\), the analytic Cauchy data
\(Y(t_*,\rho z)\) on a fixed complex disc, say \(|z|<2\), can be
made uniformly close to \(Y_0(t_*)\). This follows from their
analyticity and their value at \(z=0\).

Subtract \(Y_0(t)\), then choose common complex neighborhoods and
coefficient bounds along the compact background graph, shrinking them
if necessary.
The analytic Cauchy--Kowalevski majorant construction now gives a
positive time radius \(\delta\), independent of \(t_*\) and of the
sufficiently small \(\rho\).
More explicitly, Cauchy estimates on the fixed polydiscs give a
common coefficient majorant for the rescaled right-hand sides and
the Cauchy data. Subtracting the latter also gives a common majorant
for their first \(z\)-derivatives on a smaller disc.
The same scalar analytic Cauchy problem therefore majorizes every
one of these initial-value problems.
Thus its positive existence radius works for all of them.

At \(z=0\) the solution remains \(Y_0(t)\) by ODE uniqueness.
Restarting inside the resulting time interval, rescaling the spatial
disc again if needed, advances by a fixed positive fraction of
\(\delta\). Finitely many steps cover \([t_0,T]\).
Analytic Cauchy uniqueness identifies the overlapping germs.
Only a finite number of positive spatial-radius reductions is
required. \(\square\)

### Proposition 7.8 (nonlinear extension to every fixed finite width at the end)

For every finite \(T>0\), the solution in Proposition 7.5 admits an
analytic continuation with strictly positive curvature operator on
\[
 s\ge s_T,\qquad 0\le t\le T
 \tag{35}
\]
for some finite \(s_T\). These continuations agree on overlaps.
On each such finite-width region, its metric coefficients converge in
all derivatives to the product in (29), and \(F-3s\) converges to (30).

The quantifiers are essential: \(s_T\) may tend to infinity with \(T\).
This does not give all \(t\ge0\) at a fixed \(s\), a second smooth
singular orbit, or a complete ambient manifold.

**Proof.** For \(t>0\), introduce
\[
 u=(\log H,\log a,\log b,F-3s),\quad
 p=\partial_tu,\quad W=\xi\partial_\xi u,\qquad \xi=e^{-2s}.
\]
To make the noncharacteristic reduction explicit, write
\(u=(h,\alpha,\beta,f)\) and set
\[
 \mathcal T=h_t+\alpha_t+\beta_t+f_t,\qquad
 \Omega=\alpha_s+\beta_s-h_s+3+f_s.
\]
The four normalized diagonal equations are exactly
\[
\begin{split}
 h_{tt}&=-h_t\mathcal T-e^{-2h}\big[
 \alpha_{ss}+\beta_{ss}+f_{ss}+\alpha_s^2+\beta_s^2
       -h_s(\alpha_s+\beta_s+3+f_s)\big],\\
 \alpha_{tt}&=-\alpha_t\mathcal T
       -e^{-2h}(\alpha_{ss}+\alpha_s\Omega),\\
 \beta_{tt}&=-\beta_t\mathcal T
       -e^{-2h}(\beta_{ss}+\beta_s\Omega),\\
 f_{tt}&=-h_{tt}-\alpha_{tt}-\beta_{tt}
                  -h_t^2-\alpha_t^2-\beta_t^2.
\end{split}
\]
Here \(u_s=-2W\) and \(u_{ss}=4\xi\partial_\xi W\).
Adjoining \(u_t=p\) and \(W_t=\xi\partial_\xi p\) therefore gives
exactly a system of the form (34) for \(Y=(u,p,W)\), with analytic
coefficients. The logarithms are valid where \(H,a,b\) are nonzero.
At \(\xi=0\) the known background is the product cigar, for which
\(H=a=1\) and \(b=B(t)>0\) on every compact interval in \(t>0\).
Lemma 7.7 continues the solution from the tube to every such interval.
The relation \(W=\xi u_\xi\) persists since its difference has zero
\(t\)-derivative and vanishing initial data.
The mixed soliton residual vanishes by analytic continuation from
the overlapping initial tube. Thus this is a continuation of the
full nonlinear soliton, not merely of its linearization.

It remains to prove curvature positivity: convergence to the
semipositive product alone would not suffice.
Write the first-order coefficients as
\[
 \log H=\xi h+O(\xi^2),\quad
 \log a=-\xi y+O(\xi^2),\quad
 \log(b/B)=-\xi v+O(\xi^2),
\]
and set \(z=-h'\), \(L=B'/B\), \(D=L+\overline F'\), where
\(\overline F=2\log\cosh(t/A)\). Then \(D'=-L^2\).
The first-order equations imply
\[
 z''+Dz'=L^2z+4Lv,\qquad
 v''+Dv'=(L^2+2)v+\tfrac12Lz,
 \tag{36}
\]
with
\[
 z=2t+O(t^3),\qquad v=t^2/3+O(t^4).
\]

For clarity, (36) follows without dropping the mixed constraint.
Write \(\alpha=-y\), \(\beta=-v\), let \(\eta\) be the coefficient
of \(\xi\) in \(F-3s-\overline F\), and set
\(Q=h+\alpha+\beta+\eta\).
The linearized \(ss,tt,st,\psi\psi\) equations respectively give
\[
\begin{split}
 h''+Dh'+2h+4Q&=0,\\
 Q''+2L\beta'&=0,\\
 2Q'+h'+2L\beta&=0,\\
 \beta''+D\beta'+LQ'-2\beta&=0.
\end{split}
\]
Differentiating the first equation and substituting the third gives
the equation for \(z\); substituting the third into the fourth gives
the equation for \(v\).

Both \(L,D\) are positive on \(t>0\).
The positive integrating factor
\(\mu=\exp(\int D\,dt)\) makes (36) a cooperative system:
as long as \(z,v>0\), both \((\mu z')'\) and \((\mu v')'\) are positive.
The displayed initial expansions show \(z,z',v,v'>0\) for small
positive \(t\), so a first loss of any of these signs is impossible.
Linear ODE continuation gives these inequalities for every finite
\(t>0\).

In the five-dimensional kernel of the limiting product curvature
operator, the first-order coefficient is block diagonal.
Its \(2\times2\) block is (33), and its remaining three entries are
\[
 z'\quad\hbox{on }e_s\wedge e_t,\qquad
 4v+Lz\quad\hbox{on }e_s\wedge e_\psi,\qquad
 Ly'\quad\hbox{on }e_\phi\wedge e_\psi.
 \tag{37}
\]
These follow directly from the multiply warped curvature formulas.
All are positive for \(t>0\); their limits at zero are \(2,2,2\).
Lemma 7.6 proves positivity of the \(2\times2\) block.
The normal cigar two-plane has curvature
\(c\,\operatorname{sech}^2(t/A)>0\).
Consequently, on every fixed compact interval \(0\le t\le T\),
both this normal curvature and the least eigenvalue of the
first-order kernel block have positive lower bounds.

The analytic continuation supplies uniform Taylor remainders in
\(\xi\) on that interval. The same Schur-complement argument as in
Proposition 7.5 proves strictly positive curvature operator for all
sufficiently small \(\xi>0\), uniformly for \(0\le t\le T\).
Taking \(s_T\) large enough gives (35).
Analyticity also supplies convergence with all derivatives, and
uniqueness makes the continuations compatible. \(\square\)

Thus no obstruction occurs at a prescribed finite normal distance
when the starting point is allowed to move sufficiently far out
along the fixed surface. The unresolved region is the joint regime
in which normal distance and the axis parameter both become large,
together with the global orbit matching and completeness conditions.
This compact-time argument gives no growth rate for \(s_T\);
Proposition 7.9 uses the uniform product background to obtain a linear
bound.

### Proposition 7.9 (linear control of the continuation threshold)

For each fixed \(0<\theta<1\), there are constants \(s_*>0\) and
\(C_*>0\) such that the same analytic soliton extends with strictly
positive curvature operator throughout
\[
 s\ge s_*+C_*t,\qquad t\ge0.
 \tag{38}
\]
In particular the thresholds in Proposition 7.8 can be chosen with
\(s_T\le s_*+C_*T\) for every \(T\ge0\).

This gives a joint large-\(s,t\) region, not a complete ambient
manifold. At each fixed \(s\), (38) still allows only a finite normal
interval. It supplies no smooth matching to the second singular orbit.

**Proof.** Fix \(t_0>0\) inside the original uniform tube.
For \(t\ge t_0\), subtract the product solution \(Y_0(t)\) from the
first-order system used in Proposition 7.8 and write
\(Z=Y-Y_0\). It has the form
\[
 Z_t=\mathcal G(t,Z)+\xi\mathcal B(t,Z)Z_\xi,\qquad
 \mathcal G(t,0)=0.
 \tag{39}
\]
The coefficients are uniformly holomorphic and bounded on a fixed
complex time neighborhood of the entire half-line \(t\ge t_0\),
and for \(Z\) in a fixed small polydisc.
Indeed \(B(t)=A\tanh(t/A)\), \(B'/B\), and all required derivatives
are uniformly holomorphic away from \(t=0\) in a sufficiently narrow
strip. The logarithm of \(B\) is bounded there.
The unbounded value of the cigar potential causes no problem:
the displayed logarithmic equations depend only on its derivatives,
which are uniformly bounded. This also shows directly that subtracting
the product trajectory gives the stated uniform bounds in (39).

Here is a quantitative version of the restarting argument.
Choose a small allowed data size \(\eta>0\).
The common analytic Cauchy majorant supplies constants \(\delta>0\)
and \(M>0\) with this property: a perturbation bounded by \(\eta\)
on a fixed complex spatial disc gives a solution bounded by \(M\)
on a fixed smaller disc for a complex time interval of radius
\(\delta\). These constants work at every real \(t\ge t_0\).
Rescaling \(\xi=r z\) does not change them, because
\(\xi\partial_\xi=z\partial_z\).

All these solutions satisfy \(Z(t,0)=0\).
At a restart time inside the preceding time disc, the Schwarz bound
on the smaller spatial disc gives
\[
 |Z(t,\lambda r z)|\le C M\lambda
\]
on the next normalized data disc, where \(C\) is an absolute
disc-radius constant.
Choose a fixed \(0<\lambda<1\), smaller than the fixed Cauchy output-disc
fraction, such that \(CM\lambda\le\eta\).
Thus each restart can use exactly the fixed spatial reduction
\(r\mapsto\lambda r\), while advancing by the fixed time
step \(\delta/2\).
Induction gives a single compatible analytic continuation for all
\(t\ge t_0\). Let \(r(t)\) be the decreasing restart radius after a
fixed margin reduction. It satisfies
\[
 r(t)\ge r_*e^{-a t}
 \tag{40}
\]
for constants \(r_*,a>0\).
The bound on \(Z\) remains uniform; shrink the initial radius so that
\(r(t)\le1\).
Keeping fixed spatial and temporal margins inside each Cauchy
neighborhood, Cauchy estimates also bound its required time
derivatives and Euler derivatives \(\xi\partial_\xi\) uniformly.

Consequently the curvature matrix in the adapted orthonormal frame
is uniformly holomorphic and bounded for
\(|\xi|<c_0 r(t)\), after changing fixed constants.
The product cigar radius is bounded away from zero for \(t\ge t_0\);
the logarithmic variables keep the other metric factors bounded
away from zero on these small complex neighborhoods.
Taylor estimates therefore give, in the five-plus-one decomposition
used above,
\[
 \mathcal R(\xi,t)=
 \begin{pmatrix}
  \xi P(t)+O(\xi^2/r(t)^2)&O(\xi/r(t))\\
  O(\xi/r(t))&k_0(t)+O(\xi/r(t))
 \end{pmatrix},
 \qquad k_0(t)=c\,\operatorname{sech}^2(t/A),
 \tag{41}
\]
with uniform constants for sufficiently small \(\xi/r(t)\).

We need quantitative lower bounds for \(P(t)\), not merely its
pointwise positivity.
Set \(q=y'/y\). Lemma 7.6 gives \(0<q\le\kappa\), \(q'>0\), and
\(q'\le2\). Differentiating its Riccati equation gives
\[
 (q')'+(D+2q)q'=L^2q\ge0.
\]
Since \(D(t)\le D(t_0)\), it follows that
\(q'(t)\ge c_1e^{-K_1t}\) for positive constants.
The least eigenvalue of its \(2\times2\) curvature block is bounded
below by its determinant divided by its trace, namely
\[
 \frac{4yq'}{4+q'+q^2}\ge c_2e^{-K_1t}.
\]
We used \(y\ge2\) and the bounds on \(q,q'\).

Equation (36) also gives
\[
 z'(t)\ge c_3e^{-D(t_0)t},\qquad
 4v+Lz\ge4v(t_0)>0.
\]
Finally \(y''>0\), so \(y'(t)\ge y'(t_0)>0\), while
\(L=(2/A)\operatorname{csch}(2t/A)\ge(4/A)e^{-2t/A}\).
Thus \(Ly'\) and \(k_0\) have exponential positive lower bounds.
Combining these estimates, for some \(m,K>0\),
\[
 P(t)\ge m e^{-Kt}I,\qquad k_0(t)\ge m e^{-Kt}
 \quad(t\ge t_0).
 \tag{42}
\]

First make the last block in (41) at least
\(\tfrac12 m e^{-Kt}\).
Its Schur correction to the kernel block is then bounded above by
\(C\xi^2r(t)^{-2}e^{Kt}\).
Together with the kernel Taylor remainder, (42) shows that all
blocks are strictly positive whenever
\[
 0<\xi\le d\,r(t)^2e^{-2Kt}
 \tag{43}
\]
for a sufficiently small fixed \(d>0\).
Shrinking \(d\) also ensures the analytic-radius and last-block
conditions. By (40), it suffices to require
\[
 0<\xi\le d_*e^{-(2a+2K)t}.
\]
Since \(\xi=e^{-2s}\), this is a linear condition
\(s\ge s_*+(a+K)t\).
Increasing \(s_*\) handles the compact normal interval
\(0\le t\le t_0\) using Proposition 7.5.
All continuations agree by analytic uniqueness, so (38) describes
one exact nonlinear solution. \(\square\)

The uncontrolled part is now the complementary joint regime and
the second singular-orbit boundary, not growth of \(s_T\) faster
than every linear function. The constants in (38) are not claimed
sharp, and (38) does not imply a complete extension.

### Proposition 7.10 (a global curvature budget for a complete extension)

Suppose the analytic germ of Proposition 7.4 admits a complete steady
gradient Ricci soliton extension with the standard \(O(2)\times O(2)\) action and positive sectional
curvature. Let \(S\) be the complete reflected two-dimensional section
from Section 4, and let \(Q\) be one of its quadrants. Then
\[
 \int_Q K_S\,dA\le
 \frac{\pi}{2}-\arctan\frac{\sqrt{4/\theta-1}}3,
 \qquad
 \int_S K_S\,dA\le
 4\arctan\frac3{\sqrt{4/\theta-1}}.
 \tag{44}
\]
The second bound is strictly less than \(4\pi/3\) for \(0<\theta<1\),
and tends to zero as \(\theta\) tends to zero.

Here \(S\) is the polar section, **not** the fixed surface
\(\Sigma_p\) whose total curvature is \(2\pi\) in this construction.
Thus (44) is a necessary condition for complete extension, not the
strict fixed-surface curvature gap in (19).

**Proof.** Put \(f=F|_S\), normalized by \(f(o)=0\).
Since \(S\) is totally geodesic and \(\nabla F\) is tangent to it,
\[
 \nabla_S^2 f=\operatorname{Ric}_M|_{TS}>0.
\]
The potential is proper. Indeed, on a small compact ball about \(o\),
positive Ricci curvature gives a uniform positive lower bound for
\((F\circ\gamma)'(r_0)\) along unit-speed minimizing geodesics from
\(o\). Convexity then gives \(F(x)\ge\varepsilon d(o,x)-C\).
Properness restricts to the closed section \(S\).
The sole critical point is the nondegenerate minimum \(o\);
its regular sublevel sets are disks, and its levels are circles.
Each axis ray meets a positive level exactly once, by strict
convexity along that ray. Reflection symmetry divides the level
into four arcs.

For \(r>0\), let \(Q_r=Q\cap\{f\le r\}\), and let \(\Gamma_r\)
be its level arc. The other two boundary pieces are geodesic axis
segments. All three interior angles are \(\pi/2\): the axes are
orthogonal at \(o\), and \(\nabla f\) is tangent to them at the
level intersections. The inward geodesic curvature of \(\Gamma_r\)
is
\[
 k_r=\frac{\nabla_S^2f(T,T)}{|\nabla_S f|}>0,
\]
where \(T\) is its unit tangent. Gauss--Bonnet gives the exact identity
\[
 \int_{Q_r}K_S\,dA+\int_{\Gamma_r}k_r\,d\ell=\frac{\pi}{2}.
 \tag{45}
\]
Global minimizing properties of the axis rays are not needed.

Fix a finite \(T_0>0\). Proposition 7.8 gives, uniformly with the
needed derivatives on \(0\le t\le T_0\) as \(s\to\infty\),
\[
 h=dt^2+H(s,t)^2ds^2\longrightarrow dt^2+ds^2,\qquad
 F_s\longrightarrow3,\qquad F_t\longrightarrow u'(t),
\]
where
\[
 u(t)=2\log\cosh(t/A)-3\log2,\qquad
 d:=2/A=\sqrt{4/\theta-1}.
\]
Consequently \(F(s_r(t),t)=r\) has, for all sufficiently large \(r\),
a branch on \(0\le t\le T_0\) with \(s_r(t)\to\infty\) uniformly and
\[
 s_r'(t)\longrightarrow-u'(t)/3,\qquad
 s_r''(t)\longrightarrow-u''(t)/3.
\]
This branch is a subarc of \(\Gamma_r\), starting at the first axis.
To justify using it in a complete extension, analytic uniqueness
identifies the continuation with the normal exponential development.
The positive coefficients \(H,a,b\) on the principal part of the
finite-width end prevent local focal or singular-orbit breakdown
there. Its image remains in the open quadrant for \(t>0\).
A locally regular level arc starting at the first axis cannot wrap
around inside that quadrant: the global level there is a single
simple arc with its endpoints on the two axes. No global injectivity
of the normal coordinate system is asserted.

The curvature integral on this subarc converges to that of the
Euclidean graph \(s=\text{constant}-u(t)/3\). Thus
\[
 \lim_{r\to\infty}
 \int_{\Gamma_r[0,T_0]}k_r\,d\ell
 =
 \int_0^{T_0}
 \frac{u''(t)/3}{1+(u'(t)/3)^2}\,dt
 =
 \arctan\!\left(\frac d3\tanh(T_0/A)\right).
 \tag{46}
\]
All other portions of \(\Gamma_r\) have positive geodesic curvature.
Since \(K_S\ge0\) and the \(Q_r\) exhaust \(Q\), (45), (46), and
monotone convergence yield
\[
 \int_QK_S\,dA\le
 \frac{\pi}{2}
 -\arctan\!\left(\frac d3\tanh(T_0/A)\right).
\]
Now let \(T_0\to\infty\). Multiplying by four and using the reflections
proves (44). The limits are taken in this order; no uniform
continuation over all normal distances at a fixed \(s\) is used.
\(\square\)

This budget constrains the unknown region and second-axis matching.
It does not show that the constructed germ exceeds the budget.
No contradictory lower bound for the actual section curvature is
provided by this budget. Propositions 7.11--7.12 instead obtain an
obstruction from a change in the mixed-circle curvature sign. No
unconditional statement about the boundary topology of the
Lavoyer--Peachey family follows from (44).

### Proposition 7.11 (loss of positive curvature for the small-parameter germ)

There is \(\theta_0>0\) such that, for \(0<\theta<\theta_0\), the
specific analytic germ of Proposition 7.4 has no complete steady
gradient Ricci soliton extension with positive sectional curvature. In fact, its regular analytic
normal development already has negative mixed-circle sectional
curvature at a finite normal distance.

This excludes complete positive-curvature extension of the constructed
cigar-axis germ for small \(\theta\). It does **not** establish the
strict curvature gap for the actual Lavoyer--Peachey family, whose
axis data have not been identified with this germ.

**Proof.** Set
\[
 A=\frac2{\sqrt{4/\theta-1}},\qquad \tau=t/A,\qquad
 S(s)=\operatorname{sech}^2s,\quad w(s)=\tanh s.
\]
Fix \(s_0>0\) and a small complex \(s\)-disc around it avoiding the
zeros and poles of \(w\). All estimates below are on a smaller fixed
disc. Put
\[
 L=2\operatorname{csch}(2\tau),\quad D=2\coth(2\tau),\quad
 U=2\log\cosh\tau.
\]
The rescaled logarithmic equations in Proposition 7.8 have spatial
second derivatives multiplied by \(A^2\). After removing
\(\log A\) from \(\log b\), all coefficients and axis data are
analytic at \(A=0\). The Fuchsian majorant argument of Proposition
7.4 therefore applies with \(A\) as an additional analytic parameter,
on a common small \(\tau\)-interval. The indicial matrix is unchanged;
the rescaled free normal curvature is \(A^2N=2+A^2S/2\).
Uniqueness makes the four normalized logarithms even in \(A\).

Write
\[
 v=\bigl(\log H,\log(a/w),\log(b/(A\tanh\tau)),F-F_0-U\bigr).
\]
Its order-\(A^2\) coefficient is
\[
 v=A^2S(s)\bigl(h(\tau),h(\tau),k(\tau),f(\tau)\bigr)+O(A^4),
 \tag{47}
\]
where
\[
 h=-\tfrac12\log\cosh\tau,\qquad
 k=\tfrac14\left(\frac{\tau}{\sinh2\tau}-\tfrac12\right),\qquad
 f=\log\cosh\tau+\tfrac14\tau\tanh\tau.
 \tag{48}
\]
Here (47) is initially on a fixed small interval.
For direct verification, set \(q=2h'+k'+f'\). The coefficient
equations reduce to
\[
 h''+Dh'=-1,\qquad
 k''+Dk'+Lq=0,\qquad q'=-2Lk',
\]
and
\[
 q=\tfrac14\bigl(\coth2\tau-2\tau\operatorname{csch}^2(2\tau)\bigr).
\]
The mixed equation is \(h'+2q+2Lk=0\). These identities and the
regular axis expansions verify all the coefficient equations.
The apparent singularities at \(\tau=0\) are removable.

We next justify using (47) on a growing interval; fixed-time
parameter dependence alone would not suffice. Choose a fixed small
\(\tau_0>0\). Set \(p=v_\tau\), \(W=A v_s\), and retain the state
\(X=(v_h,p,W)\). The other three values of \(v\) do not enter the
equations, and can be recovered by integration of \(p\). The exact
first-order system has the form
\[
 X_\tau=M(\tau)X+A\,\mathcal B(\tau,s,A,X)X_s+
                    \mathcal R(\tau,s,A,X),
 \tag{49}
\]
where, uniformly for \(\tau\ge\tau_0\) on fixed complex neighborhoods,
\[
 \mathcal R=O(A^2)+O(A X)+O(X^2).
\]
This follows directly by substituting \(v_s=W/A\) into the
logarithmic equations: \(A^2v_{ss}=A W_s\), the products of two
spatial derivatives give quadratic \(W\)-terms, and cross terms
with the fixed axis derivatives are \(O(AW)\).
Also \(W_\tau=A p_s\). The coefficients are uniformly analytic:
\(D,L\) and their needed derivatives are bounded on a common narrow
complex time strip away from zero.

The \(A=0\) linear equations for \(p\) are
\[
\begin{split}
 p_h'&=-Dp_h,&p_a'&=-Dp_a,\\
 p_b'&=-(D+L)p_b-L(p_h+p_a+p_f),\\
 p_f'&=-p_h'-p_a'-p_b'-2Lp_b .
\end{split}
\]
The remaining equations are \(v_h'=p_h\), \(W'=0\).
Their forward fundamental matrices have a common bound for
\(\tau\ge\sigma\ge\tau_0\). Indeed \(p_h,p_a\) decay at least at rate \(2\). If
\(r=p_h+p_a+p_b+p_f\), the remaining subsystem is
\(p_b'=-Dp_b-Lr,\ r'=-2Lp_b\), a perturbation of
\(\operatorname{diag}(-2,0)\) by an integrable matrix.
Variation of constants and integral Gronwall give a uniform bound;
integrating the decaying \(p_h\) also bounds \(v_h\).
In particular, one must retain this bound, rather than replacing it
by an exponentially growing estimate in \(\tau\).

For completeness, the analytic estimates from (49) can be obtained
by a scalar majorant. Duhamel's formula with the bounded propagator,
followed by coefficientwise majorization in the \(s\)-disc, is bounded
by an analytic scalar transport equation
\[
 z_\tau=\frac{C\phi(x)}{(1-z/\rho)^m}
       \{A^2+Az+z^2+A(1+z)z_x\},
 \qquad z(\tau_0,x)=C A^2\phi(x),
 \tag{50}
\]
after increasing fixed constants. Here
\(\phi(x)=(1-x/r_0)^{-m}\), for a fixed integer \(m\ge1\),
majorizes the analytic coefficient bounds on a larger disc, and
\(\rho>0\) is fixed. The rational factor supplies a coefficientwise
majorant for analytic dependence on the state. It stays bounded
while \(z\) is in a sufficiently small fixed disc.
On a smaller disc, its characteristics move by \(O(A\tau)\);
along them the scalar amplitude is
\(O(A^2(1+\tau))\) for
\(\tau\le2\log(1/A)\). Both bounds follow from the scalar
differential inequalities in (50), since
\(A\log(1/A)\to0\). They keep the characteristics inside the
larger disc and the amplitude inside the chosen small polydisc.
Thus the analytic majorant constructs and bounds one compatible
solution throughout this growing interval.

Insert the approximation (47)--(48) into the retained state, including
\(W=A^3\partial_s[S(h,h,k,f)]\). Its residual in (49) is
\(O(A^4P(\tau))\) for a fixed polynomial \(P\): all order-\(A^2\)
terms cancel, the leading \(W\)-equation cancels, and (48) and their
needed derivatives have at most polynomial growth.
Its initial error at \(\tau_0\) is \(O(A^4)\).
The same Duhamel majorant for the error has source
\(O(A^4P(\tau))\), spatial derivative coefficient \(O(A)\), and
remaining linear coefficients \(O(A)+O(A^2P(\tau))\).
Characteristics and Gronwall therefore bound this error by
\(O(A^4P_1(\tau))\) on a smaller disc throughout the same interval.
Recovering the omitted logarithms by integrating \(p\), and using
Cauchy estimates for \(s\)-derivatives, only changes the polynomial.
Consequently (47) holds there with error
\(O(A^4P_2(\tau))\), including its first \(s\)- and
\(\tau\)-derivatives. The constraint \(W=A v_s\) persists by its
evolution equation. The mixed soliton equation remains zero by
analytic continuation from the initial germ.

Now take \(\tau_A=2\log(1/A)\). Equations (47)--(48) give, at \(s_0\),
\[
\begin{split}
 a_s&=S+o(1),&
 a_t&=-\tfrac12 A wS+o(A),\\
 b_s&=\tfrac14 A^3wS+o(A^3),&
 b_t&=O\!\left(A^4P_3(\log(1/A))\right).
\end{split}
\]
We used \(h'\to-1/2\), \(k\to-1/8\),
\(L(\tau_A)=O(A^4)\), and
\(k'(\tau_A)=O(A^4\log(1/A))\).
Also \(H\to1\) and \(ab/A\to w>0\). Hence
\[
 \langle\nabla a,\nabla b\rangle_h
 =\frac{a_sb_s}{H^2}+a_tb_t
 =\tfrac14A^3wS^2+o(A^3)>0.
\]
The mixed-circle curvature formula (8) yields
\[
 K(e_\phi,e_\psi)
 =-\frac{\langle\nabla a,\nabla b\rangle_h}{ab}
 =-\tfrac14 A^2\operatorname{sech}^4(s_0)+o(A^2)<0.
 \tag{51}
\]
The uniform estimates keep \(H,a,b\) nonzero along the whole
developed segment away from the initial polar axis.
In any complete steady gradient Ricci soliton extension with positive
sectional curvature, the normal geodesics
exist up to this finite distance, and analytic uniqueness identifies
their regular development with the continuation just constructed.
Thus (51) contradicts positive sectional curvature. Since \(A\to0\)
if and only if \(\theta\to0\), this proves the proposition.
\(\square\)

No explicit value of \(\theta_0\) is supplied. The following proposition treats the remaining parameters for this
germ. The strict fixed-surface gap for the actual complete examples
is not decided here.

### Proposition 7.12 (the end obstruction for every parameter)

For every \(0<\theta<1\), the specific analytic germ of Proposition
7.4 has no complete steady gradient Ricci soliton extension with
positive sectional curvature. Its regular analytic normal development
has negative mixed-circle curvature along a joint end regime
\(s,t\to\infty\).

This extends the nonextension conclusion of Proposition 7.11 to all
parameters. It remains a statement about the constructed cigar-axis
germ, not about the unknown axis data of the Lavoyer--Peachey examples.

**Proof.** Fix \(\theta\), and put
\[
 d=2/A=\sqrt{4/\theta-1}>\sqrt3,\qquad
 \kappa=\frac{\sqrt{d^2+8}-d}{2}\in(0,1),\qquad \xi=e^{-2s}.
\]
Subtract the product end and write
\[
 u=(h,\alpha,\beta,f)
   =(\log H,\log a,\log(b/B),F-3s-U),
\]
where \(B=A\tanh(t/A)\) and \(U=2\log\cosh(t/A)-3\log2\).
Let \(p=u_t\), \(W=\xi u_\xi\). At \(\xi=0\), these corrections vanish.
Choose a fixed \(t_0>0\) in the uniform tube. Its analytic end data
give convergent series in \(\xi\) at \(t_0\).

We first obtain a sharper nonlinear end estimate than the general
Cauchy-radius bound in Proposition 7.9. For the \(n\)-th linear
harmonic, set \(q=h+\alpha+\beta+f\) and \(a_n=4n^2-6n\).
The logarithmic equations give exactly
\[
\begin{split}
 \alpha''+D\alpha'+a_n\alpha&=0,\\
 \beta''+D\beta'+a_n\beta+Lq'&=0,\\
 q''+2L\beta'&=0,\\
 h''+Dh'-a_nh+4n^2q&=0,
\end{split}
 \tag{52}
\]
where \(D=d\coth(dt)\), \(L=d\operatorname{csch}(dt)\).
Introduce the rescaled parameter \(\eta=\xi e^{2t}\).
Here \(p\) retains the \(t\)-derivative at fixed \(\xi\);
the rescaled derivative of \(u\) at fixed \(\eta\) is
\(p-2\eta u_\eta\).
On coefficients in \(\eta\), the generator of (52) is shifted by
\(-2nI\). Use the mode norm
\(\|(u_n,p_n)\|_n=n|u_n|+|p_n|\).

The limiting constant shifted generators \(G_n\), obtained by
putting \(D=d,L=0\), have uniform semigroup bounds
\[
 \|e^{G_n\ell}\|_n\le
 \begin{cases}
 C e^{-\delta_1\ell},&n=1,\quad\delta_1=2-\kappa>1,\\
 C e^{-\delta_*\ell},&n\ge2,
 \end{cases}
 \qquad \ell\ge0,
 \tag{53}
\]
for some \(\delta_*>2\).
Here is a verification including uniformity in \(n\).
For \(n\ge2\), the only positive unshifted root is
\[
 \lambda_n=\frac{-d+\sqrt{d^2+4a_n}}2
\]
from the \(h\)-equation. The polynomial \(x^2+dx-a_n\) is positive
at \(x=2n-2\), since its value there is
\(2((d-1)n+2-d)>0\). Thus \(2n-\lambda_n>2\);
these gaps tend to \((d+3)/2>2\). Their infimum is strictly above
\(2\). The \(\alpha,\beta\) roots have nonpositive real parts before
the shift. Remove the forcing of \(h\) by \(q''=0\) using
\[
 \widehat h=h-c_nq-e_nq',\qquad
 c_n=4n^2/a_n,\quad e_n=dc_n/a_n.
\]
This transformation and its inverse are uniformly bounded in the
mode norm. For large \(n\), the oscillator and hyperbolic blocks
have uniformly bounded eigenvector transforms in that norm.
The remaining \(q\) block is not diagonalizable: its bound is
\(C(1+n\ell)e^{-2n\ell}\), which is uniform at any smaller fixed
decay rate. Finitely many modes, including repeated roots, are
handled separately. One may choose \(\delta_*\) strictly between
\(2\) and both \(4\) and the infimum of the displayed gaps.
For \(n=1\), the only growing roots are the simple root \(\kappa\)
in the \(\alpha,\beta\) blocks; the other blocks are stable or
polynomial before the shift. This gives the first line of (53).
Finally, the time-dependent shifted matrices differ from \(G_n\)
by an operator bounded by \(Ce^{-dt}\), uniformly in the mode norm.

We spell out the analytic estimate because a generic exponential
Cauchy bound would be insufficient. Replace each mode norm by
\[
 N_n(X)=\sup_{\ell\ge0}e^{\delta_n\ell}
                         \|e^{G_n\ell}X\|_n,
 \qquad \delta_n=\delta_1\ (n=1),\quad
        \delta_n=\delta_*\ (n\ge2).
\]
By (53), these norms are uniformly equivalent to the original mode
norms, and the constant generators are contractive at the stated
rates. The exact time-dependent linear terms add \(b(t)N_n\) to
the upper Dini derivative, with \(\int_{t_0}^\infty b(t)\,dt<\infty\).

Let
\[
 E=\sum_{n\ge1}N_n(X_n)\rho^n,\qquad
 \mathcal DE=\sum_{n\ge1}nN_n(X_n)\rho^n.
\]
In the exact system the nonlinear forcing occurs only in the
\(p\)-equations. It consists of analytic quadratic terms in
\(h,p,W\), and terms with one Euler derivative of \(W\) multiplied
by a correction vanishing at zero. Its norm is therefore bounded,
for small \(E\), by \(CE^2+CE\mathcal DE\).
There is no extra factor \(n\) from a nonlinear forcing in the
\(u\)-equation, since that equation is exactly
\(u_t=p-2\eta u_\eta\).
Consequently
\[
 D^+E\le(-\delta_1+b(t))E+CE^2+
                  \bigl(\rho'/\rho+CE\bigr)\mathcal DE.
 \tag{54}
\]
One may also absorb \(CE^2\) into \(CE\mathcal DE\), since
\(\mathcal DE\ge E\).
Initially \(E=O(\rho_0)\), and the sum over \(n\ge2\) is
\(O(\rho_0^2)\). Choose \(\rho_0\) sufficiently small and prescribe
\(\rho'/\rho=-C_1E_{\rm bd}(t)\), where
\(E_{\rm bd}(t)=C_2\rho_0 e^{-\delta_1(t-t_0)}\).
The constants absorb the finite integrating factor from \(b\).
An integral bootstrap in (54) gives \(E\le E_{\rm bd}\).
Because \(E_{\rm bd}\) is integrable, \(\rho\) has a positive limit.

These estimates can first be applied to finite coefficient
truncations. Their equations are triangular: the zero mode
vanishes and every nonlinear term is at least quadratic.
The bounds are independent of the truncation, and give convergence
of the exact analytic solution on every smaller disc for all
\(t\ge t_0\).

There is no nonlinear forcing in mode \(1\). For
\(R=\sum_{n\ge2}N_n(X_n)\rho^n\), the forcing independent of \(R\)
is therefore \(O(\rho_0^2e^{-2\delta_1(t-t_0)})\).
The same radius decrease absorbs the derivative terms containing
\(R\), giving
\[
 D^+R\le(-\delta_*+b(t)+CE)R+
                 C\rho_0^2e^{-2\delta_1(t-t_0)}.
\]
Since both \(\delta_*\) and \(2\delta_1\) exceed \(2\),
\[
 R\le C\rho_0^2e^{-2(t-t_0)}.
\]
Schwarz and Cauchy estimates on a smaller fixed \(\eta\)-disc yield
\[
 \begin{split}
 u(\xi,t)-\xi u_1(t)&=O(\xi^2e^{2t}),\\
 p(\xi,t)-\xi p_1(t)&=O(\xi^2e^{2t}),\\
 W(\xi,t)-\xi W_1(t)&=O(\xi^2e^{2t}),
 \end{split}
 \qquad |\xi|e^{2t}<r_*.
 \tag{55}
\]
The relation \(W=\xi u_\xi\) persists, and the mixed soliton
equation vanishes by analytic continuation from the initial tube.
Thus (55) controls an exact nonlinear soliton development.

Its first modes are those in Lemma 7.6 and Proposition 7.8:
\(\alpha_1=-y,\ \beta_1=-v,\ h_1'=-z\).
They satisfy
\[
 y\sim C_y e^{\kappa t},\quad y'\sim\kappa C_y e^{\kappa t},
 \qquad
 v\sim C_v e^{\kappa t},\quad v'\sim\kappa C_v e^{\kappa t},
 \qquad C_y,C_v>0,
 \tag{56}
\]
and \(h_1=O(t)\).
For clarity, the mode-\(1\) bound gives \(O(e^{\kappa t})\)
growth. In (36), the forcing of the \(z\)-equation is then
integrable, because \(\kappa<d\); hence \(z'\) is integrable and
\(z\) has a finite limit. The Riccati equation for \(y'/y\) gives
an integrable difference from \(\kappa\), establishing \(C_y>0\).
Variation of constants in the \(v\)-equation gives a finite
growing-mode coefficient, since its additional forcing and
coefficient errors have integrable growing-mode projections.
At any fixed positive time, comparison with a sufficiently small
positive multiple of \(y\) gives \(v\ge c y,\ v'\ge c y'\)
thereafter, using \(L^2v+Lz/2>0\).
Thus the coefficient \(C_v\) is positive. This proves (56).

Choose
\[
 \frac1{d+\kappa}<\gamma<\frac12,\qquad
 t=\gamma\log(1/\xi).
 \tag{57}
\]
Such a choice is possible since \(d+\kappa>2\).
Set \(\varepsilon=\xi e^{\kappa t}\) and \(r=\xi^2e^{2t}\).
Then \(\xi e^{2t}\to0\), \(r/\varepsilon\to0\), and (55) applies.
Expanding the normalized gradient product, using
\(\partial_s=-2\xi\partial_\xi\), gives
\[
 \frac{\langle\nabla a,\nabla b\rangle_h}{ab}
 =
 \xi^2(4yv+y'v')-\xi Ly'
 +O(\varepsilon r+r^2+rL)+o(\varepsilon^2).
 \tag{58}
\]
The last term includes the \(H^{-2}\) correction:
\(\xi h_1+r\to0\).
By (56), the positive term in (58) is asymptotic to
\((4+\kappa^2)C_yC_v\varepsilon^2\).
The negative term divided by \(\varepsilon^2\) tends to zero,
since it is \(O(\xi^{-1}e^{-(d+\kappa)t})\).
The error \(rL/\varepsilon^2\) is
\(O(e^{(2-d-2\kappa)t})\), which tends to zero because
\(d+2\kappa=\sqrt{d^2+8}>2\). The other errors tend to zero
relative to \(\varepsilon^2\) as well. Therefore
\[
 K(e_\phi,e_\psi)
 \sim-(4+\kappa^2)C_yC_v\,\xi^2e^{2\kappa t}<0
 \tag{59}
\]
along (57).

For each sufficiently small positive \(\xi\), the uniform estimates
keep \(H,a,b\) nonzero along the entire developed normal segment
away from its initial polar axis; the initial segment is covered
by the original tube. Thus there is no intervening focal or
orbit-collapse degeneration. In any complete steady gradient
Ricci soliton extension, analytic uniqueness identifies this
regular development with the normal exponential development of
the germ, without requiring global injectivity. Equation (59)
contradicts positive sectional curvature. \(\square\)

The constructed local solitons are therefore not counterexamples
to the complete positive-curvature problem, for any parameter.
Transferring this obstruction to the actual complete family
would require new information about its axis profile and end
expansion. The condition \(w(s)/s\to0\) alone does not supply a
bounded cigar radius, the product end, or the two specific growing
mode coefficients used above.

The local construction is not a counterexample to the complete-soliton
conjecture. It proves that strict ambient curvature, the exact soliton
PDE on a neighborhood of the complete fixed surface, and the prescribed
tip ratio still do not force the strict total-curvature gap. The
global issue is complete ambient soliton extension with the required
curvature conditions, rather than finite jets or local convergence.
Proposition 7.12 excludes such a positive-curvature steady soliton
extension for this specific germ at every parameter. The original
strict-gap question for the actual complete examples remains open. A totally geodesic three-dimensional fixed submanifold also
does not automatically inherit the steady equation: its intrinsic Ricci
tensor omits normal-curvature contributions. Three-dimensional soliton
classification cannot be applied to it on total geodesicity alone.

The weighted PDE identities do not presently give the required positive
margin, and the flat-base reduction provably removes the target examples.
The isometry-group result and the [F]-dependent ambient collapse result
do not imply the strict total-curvature bound. Equation (19) is a proved
criterion, not a proof that the criterion is satisfied.

## References

* **[LP]** L. Lavoyer and L. T. Peachey, *Steady Gradient Ricci Solitons with
  O(p) × O(q) Symmetry*, Journal of Geometric Analysis **36**, article 318
  (2026). [Published article](https://doi.org/10.1007/s12220-026-02567-x);
  [arXiv v2](https://arxiv.org/html/2509.25081v2), dated 30 July 2026.
* **[F]** H. Fang, *O(3)-symmetry of four-dimensional steady Ricci solitons*,
  [arXiv:2610.36563v1](https://arxiv.org/html/2610.36563v1), 27 October 2026,
  Theorem 1.1 and equation (1.2). Recent preprint; used only for Corollary 3.1.

* **[S]** W.-X. Shi, *Deforming the metric on complete Riemannian manifolds*,
  Journal of Differential Geometry **30** (1989), 223–301.
  [DOI](https://doi.org/10.4310/jdg/1214443292).
  The standard bounded-curvature Ricci-flow gradient estimate is used
  for the pointwise limits in Proposition 6.4.

* **[CZ]** H.-D. Cao and X.-P. Zhu, *A complete proof of the Poincaré
  and geometrization conjectures — application of the Hamilton–Perelman
  theory of the Ricci flow*, Asian Journal of Mathematics **10** (2006),
  165–492. [Text](https://www.renyi.hu/~stipsicz/perelman/cao-zhu.pdf).
  Theorems 4.1.5 and 4.2.1 supply pointed compactness and the explicit
  Gromoll–Meyer injectivity estimate; Section 2.2 gives the tensor
  strong maximum principle used in Proposition 6.6.
* **[H]** R. S. Hamilton, *Four-manifolds with positive curvature operator*,
  Journal of Differential Geometry **24** (1986), 153–179.
  [DOI](https://doi.org/10.4310/jdg/1214440433).
* **[K]** B. L. Kotschwar, *A local version of Bando's theorem on the
  real-analyticity of solutions to the Ricci flow*, Bulletin of the
  London Mathematical Society **45** (2013), 153–158.
  [DOI](https://doi.org/10.1112/blms/bds074);
  [text](https://arxiv.org/pdf/1111.0355).
  Corollary 1.3 and its following observation give analyticity of
  solitons and their vector fields.
