import os
RESULTS = os.environ.get('FTL_RESULTS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'))
os.makedirs(RESULTS, exist_ok=True)
import numpy as np, sys, json, time
import wfcheck as W

def safe_analyse(g3, h, **kw):
    try:
        return W.analyse(g3, h, h, h, **kw)
    except Exception as ex:
        print('ANALYSIS FAILED:', repr(ex)); sys.stdout.flush()
        return None
from scipy.integrate import simpson

results = []
T0 = time.time()


def total_energy_axisym(rho, xs, ys, sl):
    """Integrate an axisymmetric (about x) Eulerian density over 3D using the y>=0 half of the z=0 plane."""
    r2 = rho.copy()
    xi = xs[sl[0]]
    yi = ys[sl[1]]
    R2 = r2[sl]
    m = yi >= 0
    inner = simpson(R2[:, m] * 2 * np.pi * yi[m][None, :], x=yi[m], axis=1)
    return simpson(inner, x=xi)


def masked_stats(out):
    sl = out['slice']
    Tf = out['Tf'][(slice(None), slice(None)) + sl]
    mag = np.max(np.abs(Tf), axis=(0, 1))
    mask = mag > 1e-8 * np.nanmax(mag)
    ev = out['eigen']; smp = out['sampled']
    extra = {'region mask fraction': float(np.mean(mask)), 'max |T| [J/m^3]': float(np.nanmax(mag))}
    for k in ('NEC', 'WEC', 'SEC', 'DEC'):
        extra[k + ' min (eigen, masked) [J/m^3]'] = float(np.nanmin(np.where(mask, ev[k], np.nan)))
        extra[k + ' min (sampled, masked)'] = float(np.nanmin(np.where(mask, smp[k], np.nan)))
    extra['typeI fraction (masked)'] = float(np.mean(ev['typeI'][mask]))
    return extra


def report(label, out, xs, ys, extra=None, axisym=True):
    s = W.summary(out, label)
    s.update(masked_stats(out))
    sl = out['slice']
    rho = out['rho_E']
    rc = rho[sl]
    s['fraction of points with rho_E < 0'] = float(np.mean(rc < 0))
    i, j = np.unravel_index(np.nanargmin(rc), rc.shape)
    s['argmin rho_E (x, y) [m]'] = [float(xs[sl[0]][i]), float(ys[sl[1]][j])]
    if axisym:
        s['total Eulerian energy E = int rho_E dV [J]'] = float(total_energy_axisym(rho, xs, ys, sl))
        s['E / (M_sun c^2)'] = s['total Eulerian energy E = int rho_E dV [J]'] / (1.98892e30 * W.C ** 2)
    if extra:
        s.update(extra)
    results.append(s)
    print(json.dumps(s, indent=1)); sys.stdout.flush()
    return s


# ---------------- Alcubierre comoving ----------------
h = 0.1
xs, ys, X, Y, Z = W.grid(10.0, 10.0, h)
for v in (0.5, 0.9, 2.0):
    g3 = W.alcubierre_comoving(X, Y, Z, v, 5.0, 0.5)
    out = safe_analyse(g3, h, n_dir=400)
    if out is None:
        results.append({'metric': 'FAILED (degenerate metric, singular at some grid points)'}); continue
    # analytic total energy  E = -(v^2/12) * c^4/G * int f'^2 r^2 dr
    rr = np.linspace(0, 20, 200001)
    f = W.shape_alcubierre(rr, 5.0, 0.5)
    fp = np.gradient(f, rr)
    Ean = -(v ** 2 / 12) * W.C ** 4 / W.GN * simpson(fp ** 2 * rr ** 2, x=rr)
    report('Alcubierre comoving, v=%.1f c, R=5 m, sigma=0.5/m, h=0.1 m' % v, out, xs, ys,
           extra={'analytic total energy [J]': float(Ean)})

# ---------------- Van Den Broeck comoving ----------------
h = 0.05
xs, ys, X, Y, Z = W.grid(8.0, 8.0, h)
for as_coded in (False, True):
    g3 = W.vandenbroeck_comoving(X, Y, Z, 0.1, 2.0, 1.0, 5.0, 1.0, 0.5, as_coded=as_coded)
    out = safe_analyse(g3, h, n_dir=400)
    if out is None:
        results.append({'metric': 'FAILED (degenerate metric, singular at some grid points)'}); continue
    report('Van Den Broeck comoving, v=0.1 c, R1=2 m, s1=1/m, R2=5 m, s2=1/m, A=0.5, h=0.05 m, %s' %
           ('as coded in repo (g_tt from f)' if as_coded else 'Galilean comoving form (g_tt from 1-f)'), out, xs, ys)

# ---------------- Modified Time comoving ----------------
h = 0.1
xs, ys, X, Y, Z = W.grid(10.0, 10.0, h)
for A in (2.0, 0.5):
    for as_coded in (False, True):
        g3 = W.modified_time_comoving(X, Y, Z, 0.5, 5.0, 0.5, A, as_coded=as_coded)
        out = safe_analyse(g3, h, n_dir=400)
        if out is None:
            results.append({'metric': 'FAILED: Modified Time as coded, A=%.1f (degenerate metric, lapse^2 = 0 inside the bubble)' % A}); continue
        report('Modified Time comoving, v=0.5 c, R=5 m, sigma=0.5/m, A=%.1f, h=0.1 m, %s' %
               (A, 'as coded in repo (g_tt with f v)' if as_coded else 'Galilean comoving form (g_tt with (1-f) v)'), out, xs, ys)

# ---------------- Lentz comoving ----------------
h = 0.1
sL = 30 / 7
xs, ys, X, Y, Z = W.grid(15.0, 15.0, h)
for sm in (0, 3):
    g3 = W.lentz_comoving(X, Y, Z, 0.1, sL, smooth_cells=sm)
    out = safe_analyse(g3, h, n_dir=400)
    if out is None:
        results.append({'metric': 'FAILED (degenerate metric, singular at some grid points)'}); continue
    sl = out['slice']
    rho = out['rho_E'][sl]
    # planar (z-independent) template: energy per unit length in z over the z=0 plane
    Eperlen = simpson(simpson(rho, x=ys[sl[1]], axis=1), x=xs[sl[0]])
    report('Lentz comoving, v=0.1 c, scale=30/7 m, h=0.1 m, %s' %
           ('raw piecewise-constant template as shipped' if sm == 0 else 'template smoothed with 0.3 m Gaussian'),
           out, xs, ys, extra={'energy per unit z-length over the z=0 plane [J/m]': float(Eperlen)}, axisym=False)

# ---------------- Warp Shell comoving (Fuchs et al. 2024) ----------------
h = 0.2
R1, R2 = 10.0, 20.0
m = R2 / (2 * W.GN) * W.C ** 2 / 3
xs, ys, X, Y, Z = W.grid(30.0, 30.0, h)
worldSize = np.sqrt(2 * 30.0 ** 2 + (0.6) ** 2)
prof = W.warpshell_profiles(m, R1, R2, 0.0, 0.0, 4000, worldSize * 1.2)
print('shell mass m = %.4e kg, compactness 2GM/(R2 c^2) = %.3f' % (m, 2 * W.GN * m / (R2 * W.C ** 2)))
for v, doWarp in ((0.0, False), (0.02, True), (0.2, True), (0.5, True), (0.9, True)):
    g3 = W.warpshell_comoving(X, Y, Z, prof, v, doWarp=doWarp)
    out = safe_analyse(g3, h, n_dir=400)
    if out is None:
        results.append({'metric': 'FAILED (degenerate metric, singular at some grid points)'}); continue
    extra = {}
    report('Warp Shell comoving (Fuchs 2024), R1=10 m, R2=20 m, 2GM/R2c^2=1/3, v=%.2f c, doWarp=%d, smoothFactor=4000, h=0.2 m' % (v, doWarp),
           out, xs, ys, extra=extra)

json.dump(results, open(RESULTS + '/results.json', 'w'), indent=1)
print('done in %.0f s' % (time.time() - T0))
