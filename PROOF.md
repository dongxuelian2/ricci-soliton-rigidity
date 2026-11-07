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

The local construction is not a counterexample to the complete-soliton
conjecture. It proves that strict ambient curvature, the exact soliton
PDE on a neighborhood of the complete fixed surface, and the prescribed
tip ratio still do not force the strict total-curvature gap. The
unresolved distinction is a complete ambient extension with the required
curvature conditions, rather than finite jets or local convergence. A totally geodesic three-dimensional fixed submanifold also
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
