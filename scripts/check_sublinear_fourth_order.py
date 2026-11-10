"""Exact algebra behind Proposition 6.7; does not certify its analytic remainder.

Rebuilds the spatial source from the exact diagonal equation, including the
connection square in -H_tt/H. Uses the existing SymPy installation.
"""
import sympy as s

e,t = s.symbols("e t")
mu,ell,d = s.symbols("mu ell d", positive=True)
omega,zeta,g,j = s.symbols("omega zeta g j")
h2,h3,h4,a1,a2,a3,a4,b2,b3,b4,q2,q3,q4 = [s.Function(n)(t) for n in
    ("h2","h3","h4","a1","a2","a3","a4","b2","b3","b4","q2","q3","q4")]
h=e**2*h2+e**3*h3+e**4*h4
alpha=e*a1+e**2*a2+e**3*a3+e**4*a4
beta=e**2*b2+e**3*b3+e**4*b4
Q=e**2*q2+e**3*q3+e**4*q4
spatial=lambda f:-e**2*s.diff(f,e)
logw_prime=mu*e-omega*e**2-2*zeta*e**3
F0_prime=ell-mu*e+g*e**2+j*e**3
source=(spatial(logw_prime)+spatial(F0_prime)+spatial(spatial(Q-h))
        +(logw_prime+spatial(alpha))**2+spatial(beta)**2
        -spatial(h)*(logw_prime+F0_prime+spatial(Q-h)))
coeff=[s.expand(source).coeff(e,k) for k in range(5)]
assert s.expand(coeff[2]-mu**2)==0
assert s.expand(coeff[3]-(2*((1-mu)*omega-g)-2*mu*a1+2*ell*h2))==0
expected4=(6-4*mu)*zeta-3*j+6*q2-6*h2+(omega+a1)**2-4*mu*a2+3*ell*h3
assert s.expand(coeff[4]-expected4)==0
print("PASS: S2, S3 and S4 derived from the exact spatial source.")

log2=s.log(2)
A=-ell*mu/d
a0=2*ell*mu*log2/d**2
eta=-mu**2/d
forced_g=omega-mu**2*(1+ell**2/d**2)/(2*ell)
eta3=-2*((1-mu)*omega-g)/d
a22=-ell**2*mu/(2*d**2)
a21=(mu+ell*omega+ell*a0+ell**2*mu/d**2)/d
q22=-A**2/2
q21=ell**2*mu**2*(1+2*log2)/d**3
a20,q20,h30=s.symbols("a20 q20 h30")
poly={a1:A*t+a0, h2:eta*t+mu*a0/ell, h3:eta3*t+h30,
      a2:a22*t**2+a21*t+a20, q2:q22*t**2+q21*t+q20}
rhs4=-s.diff(h2,t)*s.diff(q2,t)-coeff[4]+2*mu**2*h2
rhs4=s.expand(rhs4.subs(poly).doit().subs(g,forced_g))
assert s.simplify(rhs4.coeff(t,2))==0
Ct=s.factor(rhs4.coeff(t,1))
assert s.simplify(Ct-mu**2*((d**2+ell**2)-mu**2*(2*d**2+ell**2))/d**3)==0
curvature=s.factor(-Ct/d-eta**2)
assert s.simplify(curvature+mu**2*(1-mu**2)*(d**2+ell**2)/d**4)==0
assert not curvature.has(omega,zeta,j,a20,q20,h30)
print("PASS: quadratic source cancellation and full fourth-order curvature coefficient:",curvature)

x,r,N2=s.symbols("x r N2",positive=True)
z=1/s.sinh(x)
y1=s.coth(x);y2=x*s.coth(x)-1
assert s.simplify(s.diff(z,x,2)-2*z**3-z)==0
assert s.simplify(y1*s.diff(y2,x)-s.diff(y1,x)*y2-1)==0
integrand=((1-r*r)/(2*r))*r*r*2/(1-r*r)
assert s.integrate(s.cancel(integrand),(r,0,1))==s.Rational(1,2)
Y=-N2*x*x/(3*d*d)
assert s.limit(Y*s.diff(z,x)-s.diff(Y,x)*z,x,0)==N2/d**2
print("PASS: adjoint identity, Wronskian, endpoint boundary term and exact integral 1/2.")
print("Uniform normal-development estimates and applicability to LP axes remain mathematical issues.")
