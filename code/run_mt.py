import numpy as np, sys, json, time
sys.path.insert(0, '/home/claude/ftl_check')
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


json.dump(results, open('/home/claude/ftl_check/results_mt.json','w'), indent=1)
