import sympy as sp
exec(open('formula_test.py').read().split('def ric_formula')[0])
subs={f:sp.Rational(3,10)*x*y+sp.Rational(1,7)*y**2, b:sp.Rational(1,3)*x**2-sp.Rational(1,5)*x*y, w:sp.Rational(1,4)*y*x}
pt={x:sp.Rational(7,10),y:sp.Rational(-3,10)}
def ric_expanded(i,j):
    e=0
    for a in range(4):
        for b_ in range(4):
            e+=sp.Rational(1,2)*ginv[a,b_]*(d2g(a,j,i,b_)+d2g(i,b_,a,j)-d2g(a,b_,i,j)-d2g(i,j,a,b_))
    for a in range(4):
        for b_ in range(4):
            for c in range(4):
                for d in range(4):
                    G2=ginv[a,c]*ginv[b_,d]
                    e+=G2*(-(Gam[c][a][d]+Gam[d][a][c])*Gam[b_][i][j] + (Gam[c][j][d]+Gam[d][j][c])*Gam[b_][i][a] + Gam[c][a][b_]*Gam[d][i][j] - Gam[c][j][b_]*Gam[d][i][a])
    return e
def ric_compact(i,j):
    e=0
    for a in range(4):
        for b_ in range(4):
            e+=sp.Rational(1,2)*ginv[a,b_]*(d2g(a,j,i,b_)+d2g(i,b_,a,j)-d2g(a,b_,i,j)-d2g(i,j,a,b_))
            for c in range(4):
                for d in range(4):
                    e+=ginv[a,b_]*ginv[c,d]*(Gam[a][i][c]*Gam[b_][j][d]-Gam[a][i][j]*Gam[b_][c][d])
    return e
for (m,n) in [(0,0),(0,1),(1,1),(1,2),(2,2),(3,3),(0,2)]:
    a1=ric_std(m,n).subs(subs).doit().subs(pt); a2=ric_expanded(m,n).subs(subs).doit().subs(pt); a3=ric_compact(m,n).subs(subs).doit().subs(pt)
    print(m,n, sp.N(a1,8), sp.N(a2,8), sp.N(a3,8))
