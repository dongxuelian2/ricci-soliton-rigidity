"""Exact identities for the relative-derivative obstruction in PROOF.md.

Asymptotic positivity and the conditional Riccati argument are proved in
the text. This checks their algebra, not an ambient soliton construction.
"""
import sympy as S

s = S.symbols('s', positive=True)
mu, ell, c, beta, eps = S.symbols('mu ell c beta eps', real=True)
delta = S.Function('delta')(s)
a = mu/s+S.diff(delta, s)/(1+delta)
k = -S.diff(a, s)-a*a
expected = mu*(1-mu)/s**2-(S.diff(delta, s, 2)
            +2*mu*S.diff(delta, s)/s)/(1+delta)
assert S.simplify(k-expected) == 0
f = ell-beta/s
U = (S.diff(f, s)-k)/2
V = (f*a-k)/2
N = (ell**2+2*c-f*f-2*S.diff(f, s)-2*f*a+2*k)/2
assert S.simplify(S.diff(f, s)-k-2*U) == 0
assert S.simplify(f*a-k-2*V) == 0
assert S.simplify(f*f+2*k+4*U+4*V+2*N-ell**2-2*c) == 0
d = eps*s**(-S.Rational(5, 2))*S.cos(s)
assert S.expand(S.diff(d, s, 3)).coeff(S.sin(s)) == (
    eps/s**S.Rational(5, 2)-S.Rational(105, 4)*eps/s**S.Rational(9, 2))
aa, b, bu = S.symbols('aa b bu')
ap = -b*aa**2
app = -bu*aa**3+2*b*b*aa**3
assert S.expand(aa**4+ap**2-aa*app-aa**4*(1-b*b+bu)) == 0
print('PASS: exact axis identities, oscillatory derivative and Riccati transformation.')
