"""Velocity bound for positive-energy warp shells (Fuchs et al. 2024 construction).

For each configuration (R1, R2, compactness C = 2GM/(R2 c^2), smoothFactor):
  1. static shell (v = 0): energy conditions, rho0, S0 (frame spatial stress)
  2. linear run (v = eps): J_i = T_(0)(i)/v, linearity check with 2*eps
  3. first-order threshold v1(x) = min_e (rho0 + S0(e,e)) / (2 J.e) over directions with J.e > 0,
     v* = min over the matter region
  4. weak-field equatorial formula v_wf = 8 pi G (rho0 + p_t0)/(c^4 |S'' + S'/r|), from the radial profile
  5. bound: v_B = 2 pi G (rho_max + p_t,max) Delta^2 / c^4, with Delta the 1% to 99% width of S
  6. bisection for the true NEC threshold v_NEC (violation when min NEC < -tol * max|T| in the matter region)
"""
import numpy as np, sys, json, time
sys.path.insert(0, '/home/claude/ftl_check')
import wfcheck as W

H = 0.2
TOL = 1e-4
NDIR = 400
EPS = 1e-3
MASK = 1e-4


def make_grid(R2, h=H):
    xmax = R2 + 6.0
    x = np.arange(-xmax, xmax + h / 2, h)
    y = np.arange(-1.0, xmax + h / 2, h)
    z = h * np.arange(-2, 3)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    return x, y, X, Y, Z


def frame_pieces(out, mask_src=None):
    sl = out['slice']
    Tf = out['Tf'][(slice(None), slice(None)) + sl]
    return Tf


def run_metric(prof, X, Y, Z, v, eigen=False):
    g3 = W.warpshell_comoving(X, Y, Z, prof, v, doWarp=(v != 0))
    return W.analyse(g3, H, H, H, n_dir=NDIR, eigen=eigen)


def matter_mask(Tf0):
    mag = np.max(np.abs(Tf0), axis=(0, 1))
    return mag > MASK * np.nanmax(mag), mag


def nec_min_ratio(out, mask, scale):
    s = out['sampled']['NEC']
    return float(np.nanmin(np.where(mask, s, np.nan)) / scale)


def analyse_config(name, R1, R2, C, smoothFactor, Rbuff=0.0, do_bisect=True, vlo=0.0, vhi=0.3):
    t0 = time.time()
    res = dict(name=name, R1=R1, R2=R2, C=C, smoothFactor=smoothFactor, Rbuff=Rbuff)
    M = C * R2 * W.C ** 2 / (2 * W.GN)
    res['M_kg'] = M
    x, y, X, Y, Z = make_grid(R2)
    worldSize = np.sqrt(2 * (R2 + 10) ** 2)
    prof = W.warpshell_profiles(M, R1, R2, Rbuff, 0.0, smoothFactor, worldSize * 1.2)

    # ---- 1. static shell
    out0 = run_metric(prof, X, Y, Z, 0.0, eigen=True)
    sl = out0['slice']
    xs, ys = x[sl[0]], y[sl[1]]
    Tf0 = frame_pieces(out0)
    mask, mag0 = matter_mask(Tf0)
    scale = float(np.nanmax(mag0))
    res['static'] = {k: float(np.nanmin(np.where(mask, out0['sampled'][k], np.nan)) / scale) for k in ('NEC', 'WEC', 'SEC', 'DEC')}
    res['static_eigen'] = {k: float(np.nanmin(np.where(mask, out0['eigen'][k], np.nan)) / scale) for k in ('NEC', 'WEC', 'SEC', 'DEC')}
    res['static_typeI_fraction'] = float(np.mean(out0['eigen']['typeI'][mask]))
    res['rho_max_J_m3'] = float(np.nanmax(Tf0[0, 0]))
    rho0 = Tf0[0, 0]
    S0 = Tf0[1:, 1:]

    # ---- 2. linear runs
    outA = run_metric(prof, X, Y, Z, EPS)
    outB = run_metric(prof, X, Y, Z, 2 * EPS)
    JA = frame_pieces(outA)[0, 1:] / EPS          # (3, N1, N2)
    JB = frame_pieces(outB)[0, 1:] / (2 * EPS)
    lin_err = float(np.nanmax(np.abs(JA - JB)[:, mask]) / np.nanmax(np.abs(JA)[:, mask]))
    res['linearity_error'] = lin_err
    J = JA
    # rho and S change at O(v^2): check
    res['rho_change_O(v2)'] = float(np.nanmax(np.abs(frame_pieces(outA)[0, 0] - rho0)[mask]) / scale)

    # ---- 3. first-order threshold
    dirs = W.fibonacci_sphere(2000)
    v1 = np.full(rho0.shape, np.inf)
    for e in dirs:
        Je = np.einsum('i...,i->...', J, e)
        num = rho0 + np.einsum('ij...,i,j->...', S0, e, e)
        with np.errstate(divide='ignore', invalid='ignore'):
            cand = np.where((Je > 0) & (num > 0), num / (2 * Je), np.inf)
        v1 = np.minimum(v1, cand)
    v1m = np.where(mask, v1, np.inf)
    i, j = np.unravel_index(np.argmin(v1m), v1m.shape)
    res['v_first_order'] = float(v1m[i, j])
    res['v_first_order_location'] = dict(x=float(xs[i]), y=float(ys[j]), r=float(np.hypot(xs[i], ys[j])))
    # also the simple estimate using e parallel to J
    Jn = np.sqrt(np.einsum('i...,i...->...', J, J))
    with np.errstate(divide='ignore', invalid='ignore'):
        Jhat = J / Jn
        num_par = rho0 + np.einsum('ij...,i...,j...->...', S0, Jhat, Jhat)
        v1_par = np.where(num_par > 0, num_par / (2 * Jn), np.inf)
    res['v_first_order_e_parallel_J'] = float(np.nanmin(np.where(mask & (Jn > 0), v1_par, np.inf)))

    # ---- 4. weak-field equatorial formula from the radial profile
    r = prof['r']; S = prof['S']
    dr = r[1] - r[0]
    Sp = np.gradient(S, dr); Spp = np.gradient(Sp, dr)
    # rho0 and tangential pressure along the +y axis (x = 0 column), as functions of r = y
    ix0 = np.argmin(np.abs(xs - 0.0))
    rho_line = rho0[ix0, :]; pt_line = S0[0, 0, ix0, :]   # at x=0 on the y axis, x-direction is tangential
    yy = ys
    sel = (yy > R1 - 2) & (yy < R2 + 2)
    Spp_y = np.interp(yy, r, Spp); Sp_y = np.interp(yy, r, Sp)
    with np.errstate(divide='ignore', invalid='ignore'):
        denom = np.abs(Spp_y + Sp_y / np.where(yy > 0, yy, 1.0))
        v_wf_line = 8 * np.pi * W.GN * (rho_line + pt_line) / W.C ** 4 / denom
    good = sel & (denom > 1e-3 * np.nanmax(denom)) & mask[ix0, :]
    res['v_weak_field_equator'] = float(np.nanmin(v_wf_line[good]))
    # also weak-field J magnitude vs numerical J at the equator
    Jwf = (np.abs(Spp_y + Sp_y / np.where(yy > 0, yy, 1.0)) / (16 * np.pi)) * W.C ** 4 / W.GN
    Jnum_line = Jn[ix0, :]
    gg = good & (Jwf > 0.3 * np.nanmax(Jwf))
    res['J_numerical_over_weak_field_equator_median'] = float(np.nanmedian(Jnum_line[gg] / Jwf[gg]))

    # ---- 5. bound
    inside = (S > 0.01) & (S < 0.99)
    Delta = float(r[inside].max() - r[inside].min())
    res['Delta_m'] = Delta
    res['max_abs_Spp'] = float(np.nanmax(np.abs(Spp)))
    res['four_over_Delta2'] = 4 / Delta ** 2
    rho_max = float(np.nanmax(rho0[mask]))
    pt_max = float(np.nanmax(S0[0, 0][mask]))
    res['pt_max_over_rho_max'] = pt_max / rho_max
    res['v_bound'] = float(2 * np.pi * W.GN * (rho_max + max(pt_max, 0.0)) * Delta ** 2 / W.C ** 4)
    res['v_bound_uniform_shell'] = float(0.75 * C * (1 - R1 / R2) * (1 + max(pt_max, 0.0) / rho_max))

    # ---- 6. bisection for the true NEC threshold
    if do_bisect:
        def violates(v):
            o = run_metric(prof, X, Y, Z, v)
            return nec_min_ratio(o, mask, scale) < min(res['static']['NEC'], 0.0) - TOL
        lo, hi = vlo, vhi
        # make sure hi violates
        while not violates(hi) and hi < 0.95:
            lo, hi = hi, min(hi * 2, 0.95)
        for _ in range(8):
            mid = 0.5 * (lo + hi)
            if violates(mid):
                hi = mid
            else:
                lo = mid
        res['v_NEC_bisection'] = [lo, hi]
    res['runtime_s'] = time.time() - t0
    print(json.dumps(res, indent=1)); sys.stdout.flush()
    return res


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    configs = [
        ('base: C=1/3, R1/R2=0.5, sF=4000', 10.0, 20.0, 1 / 3, 4000),
        ('C=0.10', 10.0, 20.0, 0.10, 4000),
        ('C=0.50', 10.0, 20.0, 0.50, 4000),
        ('C=0.70', 10.0, 20.0, 0.70, 4000),
        ('R1/R2=0.25', 5.0, 20.0, 1 / 3, 4000),
        ('R1/R2=0.75', 15.0, 20.0, 1 / 3, 4000),
        ('sF=2000 (sharper shift)', 10.0, 20.0, 1 / 3, 2000),
        ('sF=8000 (smoother shift)', 10.0, 20.0, 1 / 3, 8000),
    ]
    if which != 'all':
        configs = [c for c in configs if which in c[0]]
    allres = []
    for c in configs:
        allres.append(analyse_config(*c))
        json.dump(allres, open('/home/claude/ftl_check/vbound_results.json', 'w'), indent=1)
    print('ALL DONE')
