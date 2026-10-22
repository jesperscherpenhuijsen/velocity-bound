import numpy as np, sympy as sp, time, sys
import wfcheck as W

# ---------- 1. Exact Einstein tensor for Alcubierre comoving via Christoffels (sympy) ----------
t0 = time.time()
x, y, z = sp.symbols('x y z', real=True)
v, R, s = sp.symbols('v R sigma', positive=True)
r = sp.sqrt(x**2 + y**2 + z**2)
f = (sp.tanh(s*(R + r)) + sp.tanh(s*(R - r))) / (2*sp.tanh(R*s))
beta = v*(1 - f)
X = [sp.Symbol('T'), x, y, z]           # x^0 = c t
g = sp.Matrix([[-1 + beta**2, beta, 0, 0], [beta, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
ginv = g.inv()
# Christoffel symbols Gamma^a_{bc}
Gam = [[[sum(ginv[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(4))/2
         for c in range(4)] for b in range(4)] for a in range(4)]
def ricci(b, d):
    expr = 0
    for a in range(4):
        expr += sp.diff(Gam[a][b][d], X[a]) - sp.diff(Gam[a][b][a], X[d])
        for e in range(4):
            expr += Gam[a][a][e]*Gam[e][b][d] - Gam[a][d][e]*Gam[e][b][a]
    return expr
Ric = sp.Matrix(4, 4, lambda i, j: ricci(i, j) if i <= j else 0)
for i in range(4):
    for j in range(i):
        Ric[i, j] = Ric[j, i]
Rs = sum(ginv[i, j]*Ric[i, j] for i in range(4) for j in range(4))
G = Ric - Rs*g/2
print('symbolic done in %.1f s' % (time.time()-t0)); sys.stdout.flush()
Gfun = sp.lambdify((x, y, z, v, R, s), [[G[i, j] for j in range(4)] for i in range(4)], 'numpy', cse=True)

# ---------- 2. Numerical pipeline on the same metric ----------
vv, RR, ss, h = 0.5, 5.0, 0.5, 0.1
xs, ys, Xg, Yg, Zg = W.grid(9.0, 9.0, h)
g3 = W.alcubierre_comoving(Xg, Yg, Zg, vv, RR, ss)
out = W.analyse(g3, h, h, h, n_dir=100, eigen=False)
sl = out['slice']
Gnum = out['T'] / W.KAPPA
Xc, Yc = Xg[:, :, 2], Yg[:, :, 2]
Gex = np.array(Gfun(Xc, Yc, 0*Xc, vv, RR, ss), dtype=float)
Gex = np.broadcast_to(Gex, Gnum.shape) if Gex.shape != Gnum.shape else Gex
err = np.nanmax(np.abs(Gnum[(slice(None), slice(None))+sl] - Gex[(slice(None), slice(None))+sl]))
scale = np.nanmax(np.abs(Gex[(slice(None), slice(None))+sl]))
print('Alcubierre: max|G_num - G_exact| = %.3e, max|G_exact| = %.3e, relative %.2e' % (err, scale, err/scale))

# Eulerian energy density vs Alcubierre's closed form rho = -v^2 (y^2+z^2) f'^2 / (32 pi r^2) (G=c=1)
rr = np.sqrt(Xc**2 + Yc**2)
fp = sp.lambdify((x, y, z, v, R, s), sp.diff(f, x)*x/r + sp.diff(f, y)*y/r + sp.diff(f, z)*z/r, 'numpy')  # df/dr via chain
dfdr = fp(Xc, Yc, 0*Xc, vv, RR, ss)
rho_exact = -vv**2 * Yc**2 * dfdr**2 / (32*np.pi*rr**2) * W.C**4/W.GN
rho_num = out['rho_E']
e2 = np.nanmax(np.abs(rho_num[sl] - rho_exact[sl])) / np.nanmax(np.abs(rho_exact[sl]))
print('Alcubierre: rho_E relative error vs closed form: %.2e ; rho_E min = %.3e J/m^3, max = %.3e J/m^3' % (e2, np.nanmin(rho_num[sl]), np.nanmax(rho_num[sl])))

# Hamiltonian-constraint cross-check for unit lapse, flat slices: 16 pi rho = K^2 - K_ij K^ij, K_ij = (d_i b_j + d_j b_i)/2
bx = vv*(1 - np.array(sp.lambdify((x, y, z, v, R, s), f, 'numpy')(Xg, Yg, Zg, vv, RR, ss)))
def d(F, ax):
    return np.gradient(F, h, axis=ax, edge_order=2)
dbx = [d(bx, 0)[:, :, 2], d(bx, 1)[:, :, 2], d(bx, 2)[:, :, 2]]
K = np.zeros((3, 3) + Xc.shape)
for i in range(3):
    K[0, i] += dbx[i]/2; K[i, 0] += dbx[i]/2
theta = np.trace(K)
rho_ham = (theta**2 - np.einsum("ij...,ij...->...", K, K)) / (16*np.pi) * W.C**4/W.GN
e3 = np.nanmax(np.abs(rho_num[sl] - rho_ham[sl])) / np.nanmax(np.abs(rho_ham[sl]))
print('Alcubierre: rho_E vs Hamiltonian constraint (2nd-order gradient): relative %.2e' % e3)

# ---------- 3. Schwarzschild, isotropic coordinates: vacuum test ----------
rs_ = 0.2
xs2, ys2, X2, Y2, Z2 = W.grid(6.0, 6.0, 0.1, z0=0.0)
rho_iso = np.sqrt(X2**2 + Y2**2 + Z2**2)
rho_iso = np.where(rho_iso < 1e-9, 1e-9, rho_iso)
psi = 1 + rs_/(4*rho_iso)
gtt = -((1 - rs_/(4*rho_iso))/(1 + rs_/(4*rho_iso)))**2
g3s = {(0, 0): gtt, (1, 1): psi**4, (2, 2): psi**4, (3, 3): psi**4}
outs = W.analyse(g3s, 0.1, 0.1, 0.1, n_dir=50, eigen=False)
Gs = outs['T'] / W.KAPPA
Rsq = np.sqrt(X2[:, :, 2]**2 + Y2[:, :, 2]**2)
mask = (Rsq > 1.0)
mask[:4, :] = False; mask[-4:, :] = False; mask[:, :4] = False; mask[:, -4:] = False
# compare against the size of the second derivatives (scale ~ rs/r^3)
print('Schwarzschild (r>1, rs=0.2): max|G| = %.3e, typical curvature scale rs/r^3 at r=1: %.3e' % (np.nanmax(np.abs(Gs[:, :, mask])), rs_))
print('Minkowski check: max|G| = %.3e' % np.nanmax(np.abs(W.analyse(W.minkowski_slab(X2), .1, .1, .1, n_dir=10, eigen=False)['T'])))
