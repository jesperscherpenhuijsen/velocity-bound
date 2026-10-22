"""Exact linearized momentum density J for beta_i = -v S(r) delta_i^x on the static background
   ds^2 = -e^{2a} dt^2 + e^{2b} dr^2 + r^2 dOmega^2  (x-axis = polar axis).
   8 pi J_i = D_j (Kt^j_i - delta^j_i Kt),  Kt_ij = (D_i bt_j + D_j bt_i) / (2 e^a),  bt = -S(r) dx.
"""
import sympy as sp, numpy as np, json, sys
sys.path.insert(0, '/home/claude/ftl_check')

r, th, ph = sp.symbols('r theta phi', positive=True)
a = sp.Function('a')(r); b = sp.Function('b')(r); S = sp.Function('S')(r)
coords = [r, th, ph]
gam = sp.diag(sp.exp(2*b), r**2, r**2*sp.sin(th)**2)
gaminv = gam.inv()
# one-form bt = -S dx, dx = cos(th) dr - r sin(th) dth
bt = [-S*sp.cos(th), S*r*sp.sin(th), 0]
# Christoffels of gamma
Gam = [[[sum(gaminv[k, l]*(sp.diff(gam[l, i], coords[j]) + sp.diff(gam[l, j], coords[i]) - sp.diff(gam[i, j], coords[l])) for l in range(3))/2
         for j in range(3)] for i in range(3)] for k in range(3)]
def Dcov(i, j):  # D_i bt_j
    return sp.diff(bt[j], coords[i]) - sum(Gam[k][i][j]*bt[k] for k in range(3))
Kt = sp.Matrix(3, 3, lambda i, j: (Dcov(i, j) + Dcov(j, i))/(2*sp.exp(a)))
Ktmix = gaminv*Kt            # Kt^j_i  (row j, col i)
trK = sum(Ktmix[i, i] for i in range(3))
Tmix = Ktmix - trK*sp.eye(3) # T^j_i
# divergence D_j T^j_i = d_j T^j_i + Gam^j_jk T^k_i - Gam^k_ji T^j_k
def divT(i):
    e = 0
    for j in range(3):
        e += sp.diff(Tmix[j, i], coords[j])
        for k in range(3):
            e += Gam[j][j][k]*Tmix[k, i] - Gam[k][j][i]*Tmix[j, k]
    return e
J_r = sp.simplify(divT(0)/(8*sp.pi))
J_th = sp.simplify(divT(1)/(8*sp.pi))
J_ph = sp.simplify(divT(2)/(8*sp.pi))
# orthonormal components
Jrh = sp.simplify(J_r*sp.exp(-b)); Jth_h = sp.simplify(J_th/r)
print('J_phi =', J_ph)
print('J_rhat =', Jrh)
print('J_thetahat =', Jth_h)
# weak-field check: a = b = 0
wf_r = sp.simplify(Jrh.subs({a: 0, b: 0}).doit()); wf_t = sp.simplify(Jth_h.subs({a: 0, b: 0}).doit())
print('weak field J_rhat =', wf_r)
print('weak field J_thetahat =', wf_t)
# expected weak field from Cartesian formula: 8 pi J_x = -1/2 [S'' sin^2 + S'(1+cos^2)/r], 8pi J_y = 1/2 nx ny (S'' - S'/r)
# radial/polar projections in the x-z (meridional) plane: J_rhat = J_x cos + J_perp sin, J_thetahat = -J_x sin + J_perp cos with J_perp the component along the transverse unit vector
Sp, Spp = sp.diff(S, r), sp.diff(S, r, 2)
Jx = -(Spp*sp.sin(th)**2 + Sp/r*(1+sp.cos(th)**2))/(16*sp.pi)
Jperp = sp.cos(th)*sp.sin(th)*(Spp - Sp/r)/(16*sp.pi)     # n_x n_y (S''-S'/r)/2/(8pi) with n_y = sin(th)
exp_r = sp.simplify(Jx*sp.cos(th) + Jperp*sp.sin(th)); exp_t = sp.simplify(-Jx*sp.sin(th) + Jperp*sp.cos(th))
print('weak-field radial agrees with Cartesian formula:', sp.simplify(wf_r - exp_r) == 0)
print('weak-field polar  agrees with Cartesian formula:', sp.simplify(wf_t - exp_t) == 0)

# ---- lambdify with symbols for the profile and its derivatives
S0, S1, S2, a0, a1, b0, b1 = sp.symbols('S0 S1 S2 a0 a1 b0 b1')
rep = {sp.Derivative(S, (r, 2)): S2, sp.Derivative(S, r): S1, sp.Derivative(a, r): a1, sp.Derivative(b, r): b1}
def prep(e):
    e = e.doit()
    for k, v in rep.items():
        e = e.subs(k, v)
    return e.subs({S: S0, a: a0, b: b0})
fJr = sp.lambdify((r, th, S0, S1, S2, a0, a1, b0, b1), prep(Jrh), 'numpy')
fJt = sp.lambdify((r, th, S0, S1, S2, a0, a1, b0, b1), prep(Jth_h), 'numpy')
sp.init_printing()
with open('/home/claude/ftl_check/exactJ_expressions.txt', 'w') as f:
    f.write('J_rhat = ' + sp.latex(prep(Jrh)) + '\n\nJ_thetahat = ' + sp.latex(prep(Jth_h)) + '\n')
    f.write('\nequator (theta=pi/2): J_thetahat = ' + sp.latex(sp.simplify(prep(Jth_h).subs(th, sp.pi/2))) + '\n')
    f.write('equator: J_rhat = ' + sp.latex(sp.simplify(prep(Jrh).subs(th, sp.pi/2))) + '\n')
print('equator J_thetahat =', sp.simplify(prep(Jth_h).subs(th, sp.pi/2)))
print('equator J_rhat     =', sp.simplify(prep(Jrh).subs(th, sp.pi/2)))

# ---- numerical evaluation on the shell profiles
import wfcheck as W
def shell_threshold_1d(R1, R2, C, sF):
    M = C*R2*W.C**2/(2*W.GN)
    prof = W.warpshell_profiles(M, R1, R2, 0.0, 0.0, sF, 1.2*np.sqrt(2*(R2+10)**2))
    rr = prof['r']; Sv = prof['S']; A = prof['A']; B = prof['B']
    sel = (rr > 0.5)
    rr, Sv, A, B = rr[sel], Sv[sel], A[sel], B[sel]
    av = 0.5*np.log(-A); bv = 0.5*np.log(B)
    dr = rr[1]-rr[0]
    S1v = np.gradient(Sv, dr); S2v = np.gradient(S1v, dr)
    a1v = np.gradient(av, dr); b1v = np.gradient(bv, dr); a2v = np.gradient(a1v, dr)
    # geometric-unit static shell quantities (1/m^2), then SI
    e2b = np.exp(-2*bv)
    rho = (1 - e2b*(1 - 2*rr*b1v))/(8*np.pi*rr**2)
    pt = e2b*(a2v + a1v**2 - a1v*b1v + (a1v - b1v)/rr)/(8*np.pi)
    pr = (e2b*(2*a1v/rr + 1/rr**2) - 1/rr**2)/(8*np.pi)
    # J at equator and on the axis (orthonormal, 1/m^2)
    Jt_eq = fJt(rr, np.pi/2, Sv, S1v, S2v, av, a1v, bv, b1v)
    Jr_eq = fJr(rr, np.pi/2, Sv, S1v, S2v, av, a1v, bv, b1v)
    Jr_ax = fJr(rr, 1e-9, Sv, S1v, S2v, av, a1v, bv, b1v)
    Jeq = np.sqrt(Jt_eq**2 + Jr_eq**2)
    wf_eq = np.abs(S2v + S1v/rr)/(16*np.pi)
    wall = (rr > R1 - 1) & (rr < R2 + 1) & (Jeq > 1e-3*np.nanmax(Jeq))
    with np.errstate(divide='ignore', invalid='ignore'):
        v_eq = (rho + pt)/(2*Jeq)          # equator, e along theta-hat (tangential)
        v_ax = (rho + pr)/(2*np.abs(Jr_ax)) # axis, e radial
        v_wf = (rho + pt)/(2*wf_eq)
    i = np.nanargmin(np.where(wall & (rho + pt > 0), v_eq, np.inf))
    out = dict(C=C, R1=R1, R2=R2, sF=sF,
               v_exact_equator=float(v_eq[i]), r_min=float(rr[i]),
               v_weakfield_equator=float(np.nanmin(np.where(wall & (rho+pt > 0), v_wf, np.inf))),
               v_axis_min=float(np.nanmin(np.where(wall & (rho+pr > 0) & (np.abs(Jr_ax) > 1e-3*np.nanmax(np.abs(Jr_ax))), v_ax, np.inf))),
               J_exact_over_wf_at_min=float(Jeq[i]/wf_eq[i]),
               J_exact_over_wf_median_wall=float(np.nanmedian((Jeq/wf_eq)[wall & (wf_eq > 0.3*np.nanmax(wf_eq))])))
    return out

res = []
for (R1, R2, C, sF) in [(10, 20, 1/3, 4000), (10, 20, 0.10, 4000), (10, 20, 0.50, 4000), (10, 20, 0.70, 4000),
                        (5, 20, 1/3, 4000), (15, 20, 1/3, 4000), (10, 20, 1/3, 2000), (10, 20, 1/3, 8000)]:
    o = shell_threshold_1d(R1, R2, C, sF); res.append(o)
    print('C=%.2f R1/R2=%.2f sF=%d: v1D_exact=%.4f at r=%.1f | v_wf=%.4f | axis min %.4f | J/Jwf at min %.3f, median %.3f' % (
        C, R1/R2, sF, o['v_exact_equator'], o['r_min'], o['v_weakfield_equator'], o['v_axis_min'], o['J_exact_over_wf_at_min'], o['J_exact_over_wf_median_wall']))
json.dump(res, open('/home/claude/ftl_check/exactJ_results.json', 'w'), indent=1)
