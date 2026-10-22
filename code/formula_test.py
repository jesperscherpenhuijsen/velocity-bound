import sympy as sp, itertools
x,y,z,T=sp.symbols('x y z T',real=True)
X=[T,x,y,z]
# generic-ish metric: Alcubierre-like plus a conformal factor, all functions of x,y
f=sp.Function('f')(x,y); b=sp.Function('b')(x,y); w=sp.Function('w')(x,y)
g=sp.Matrix([[-1+b**2, b, w, 0],[b, 1+f, 0, 0],[w,0,1,0],[0,0,0,1+f]])
ginv=g.inv()
def dg(a,m,n): return sp.diff(g[m,n],X[a])
def d2g(a,c,m,n): return sp.diff(g[m,n],X[a],X[c])
Gam=[[[sp.Rational(1,2)*(dg(m,a,c)+dg(c,a,m)-dg(a,m,c)) for c in range(4)] for m in range(4)] for a in range(4)]  # Gam[a][m][c] lowered
GamU=[[[sum(ginv[a,d]*Gam[d][m][c] for d in range(4)) for c in range(4)] for m in range(4)] for a in range(4)]
def ric_std(b_,d_):
    e=0
    for a in range(4):
        e+=sp.diff(GamU[a][b_][d_],X[a])-sp.diff(GamU[a][b_][a],X[d_])
        for q in range(4):
            e+=GamU[a][a][q]*GamU[q][b_][d_]-GamU[a][d_][q]*GamU[q][b_][a]
    return e
def ric_formula(m,n):
    e=0
    for a in range(4):
        for b_ in range(4):
            e+=sp.Rational(1,2)*ginv[a,b_]*(d2g(a,n,m,b_)+d2g(m,b_,a,n)-d2g(a,b_,m,n)-d2g(m,n,a,b_))
            for c in range(4):
                for d in range(4):
                    e+=ginv[a,b_]*ginv[c,d]*(Gam[a][m][c]*Gam[b_][n][d]-Gam[m][n][c]*Gam[a][b_][d])
    return e
subs={f:sp.Rational(3,10)*x*y+sp.Rational(1,7)*y**2, b:sp.Rational(1,3)*x**2-sp.Rational(1,5)*x*y, w:sp.Rational(1,4)*y*x}
pt={x:sp.Rational(7,10),y:sp.Rational(-3,10)}
for (m,n) in [(0,0),(0,1),(1,1),(1,2),(2,2),(3,3),(0,2)]:
    a1=ric_std(m,n).subs(subs).doit().subs(pt); a2=ric_formula(m,n).subs(subs).doit().subs(pt)
    print(m,n, sp.N(a1,8), sp.N(a2,8))
