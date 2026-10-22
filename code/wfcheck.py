"""Independent re-implementation of the Warp Factory pipeline (Helmerich et al. 2024)
for stationary (comoving) metrics on a thin z-slab.

Conventions: signature (-,+,+,+); coordinates x^0 = c t, x^1..3 = x, y, z in metres.
Metric components are dimensionless; T_{mu nu} = c^4/(8 pi G) G_{mu nu} in J/m^3.
All metrics are evaluated on 5 z-planes; derivatives are 4th-order central differences
and all curvature quantities are returned on the central plane z = z0 only.
"""
import numpy as np

C = 2.99792458e8
GN = 6.67430e-11
KAPPA = C**4 / (8 * np.pi * GN)  # J/m^3 per unit Einstein tensor component


# ----------------------------------------------------------------------------
# finite differences on a (Nx, Ny, 5) array, result on the central z-plane
# ----------------------------------------------------------------------------
def d1(F, axis, h):
    """4th-order central first derivative along axis 0 or 1 of a 2D array (interior only)."""
    out = np.full_like(F, np.nan)
    if axis == 0:
        out[2:-2, :] = (-(F[4:, :] - F[:-4, :]) + 8 * (F[3:-1, :] - F[1:-3, :])) / (12 * h)
    else:
        out[:, 2:-2] = (-(F[:, 4:] - F[:, :-4]) + 8 * (F[:, 3:-1] - F[:, 1:-3])) / (12 * h)
    return out


def d2(F, axis, h):
    out = np.full_like(F, np.nan)
    if axis == 0:
        out[2:-2, :] = (-F[4:, :] + 16 * F[3:-1, :] - 30 * F[2:-2, :] + 16 * F[1:-3, :] - F[:-4, :]) / (12 * h * h)
    else:
        out[:, 2:-2] = (-F[:, 4:] + 16 * F[:, 3:-1] - 30 * F[:, 2:-2] + 16 * F[:, 1:-3] - F[:, :-4]) / (12 * h * h)
    return out


def slab_derivs(g3, hx, hy, hz):
    """g3: dict (mu,nu)->array (Nx,Ny,5). Returns g, dg, d2g on the central plane.
    Index 0 = time (all time derivatives zero: stationary metric)."""
    k = 2
    N1, N2 = g3[(0, 0)].shape[:2]
    g = np.zeros((4, 4, N1, N2))
    dg = np.zeros((4, 4, 4, N1, N2))      # dg[a, m, n] = d_a g_{mn}
    d2g = np.zeros((4, 4, 4, 4, N1, N2))  # d2g[a, b, m, n] = d_a d_b g_{mn}
    for (m, n), A in g3.items():
        gc = A[:, :, k]
        gz1 = (-(A[:, :, k + 2] - A[:, :, k - 2]) + 8 * (A[:, :, k + 1] - A[:, :, k - 1])) / (12 * hz)
        gz2 = (-A[:, :, k + 2] + 16 * A[:, :, k + 1] - 30 * A[:, :, k] + 16 * A[:, :, k - 1] - A[:, :, k - 2]) / (12 * hz * hz)
        gx = d1(gc, 0, hx)
        gy = d1(gc, 1, hy)
        gxx = d2(gc, 0, hx)
        gyy = d2(gc, 1, hy)
        gxy = d1(gx, 1, hy)
        gxz = d1(gz1, 0, hx)
        gyz = d1(gz1, 1, hy)
        for (p, q) in ((m, n), (n, m)):
            g[p, q] = gc
            dg[1, p, q], dg[2, p, q], dg[3, p, q] = gx, gy, gz1
            d2g[1, 1, p, q], d2g[2, 2, p, q], d2g[3, 3, p, q] = gxx, gyy, gz2
            d2g[1, 2, p, q] = d2g[2, 1, p, q] = gxy
            d2g[1, 3, p, q] = d2g[3, 1, p, q] = gxz
            d2g[2, 3, p, q] = d2g[3, 2, p, q] = gyz
    return g, dg, d2g


# ----------------------------------------------------------------------------
# curvature
# ----------------------------------------------------------------------------
def einstein_from_derivs(g, dg, d2g):
    """Ricci tensor from first and second metric derivatives (standard coordinate formula),
    then Einstein tensor. Returns ginv, R_{mn}, R, G_{mn}."""
    N1, N2 = g.shape[2:]
    gm = np.moveaxis(g, (0, 1), (2, 3))              # (N1,N2,4,4)
    ginv = np.moveaxis(np.linalg.inv(gm), (2, 3), (0, 1))
    # Christoffel of the first kind: Gam[a,m,c] = 1/2 (d_m g_{ac} + d_c g_{am} - d_a g_{mc})
    Gam = 0.5 * (np.einsum('mac...->amc...', dg) + np.einsum('cam...->amc...', dg) - np.einsum('amc...->amc...', dg))
    # second-derivative part
    t1 = np.einsum('ab...,anmb...->mn...', ginv, d2g)   # g^{ab} d_a d_n g_{mb}
    t2 = np.einsum('ab...,mban...->mn...', ginv, d2g)   # g^{ab} d_m d_b g_{an}
    t3 = np.einsum('ab...,abmn...->mn...', ginv, d2g)   # g^{ab} d_a d_b g_{mn}
    t4 = np.einsum('ab...,mnab...->mn...', ginv, d2g)   # g^{ab} d_m d_n g_{ab}
    Ric = 0.5 * (t1 + t2 - t3 - t4)
    # Gamma-Gamma part (verified symbolically against the Christoffel definition, formula_test2.py):
    # g^{ab} g^{cd} (Gam_{amc} Gam_{bnd} - Gam_{amn} Gam_{bcd}), Gam with first index lowered
    p1 = np.einsum('amc...,ab...,bnd...,cd...->mn...', Gam, ginv, Gam, ginv, optimize=True)
    p2 = np.einsum('amn...,ab...,bcd...,cd...->mn...', Gam, ginv, Gam, ginv, optimize=True)
    Ric = Ric + p1 - p2
    Rs = np.einsum('mn...,mn...->...', ginv, Ric)
    Gmn = Ric - 0.5 * Rs * g
    return ginv, Ric, Rs, Gmn


# ----------------------------------------------------------------------------
# Eulerian orthonormal frame and energy conditions
# ----------------------------------------------------------------------------
def eulerian_frame(g):
    """Return tetrad e[a, mu] (a = frame index) for the Eulerian observer of the t = const slices."""
    N1, N2 = g.shape[2:]
    gam = np.moveaxis(g[1:, 1:], (0, 1), (2, 3))          # (N1,N2,3,3)
    L = np.linalg.cholesky(gam)
    Linv = np.linalg.inv(L)
    gaminv = np.linalg.inv(gam)
    beta_lo = np.moveaxis(g[0, 1:], 0, 2)                  # (N1,N2,3)
    beta_up = np.einsum('...ij,...j->...i', gaminv, beta_lo)
    alpha2 = np.einsum('...i,...i->...', beta_up, beta_lo) - g[0, 0]
    alpha = np.sqrt(alpha2)
    e = np.zeros((N1, N2, 4, 4))
    e[..., 0, 0] = 1.0 / alpha
    e[..., 0, 1:] = -beta_up / alpha[..., None]
    e[..., 1:, 1:] = Linv                                 # e_(a)^i = (L^{-1})_{a i}
    return np.moveaxis(e, (2, 3), (0, 1)), alpha, beta_up  # e[a, mu, N1, N2]


def frame_tensor(Tmn, e):
    """T_(a)(b) = T_{mu nu} e_(a)^mu e_(b)^nu"""
    return np.einsum('mn...,am...,bn...->ab...', Tmn, e, e)


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)


def energy_conditions_sampled(Tf, n_dir=400, etas=np.linspace(0, 2.0, 10)):
    """Minimum over sampled null/timelike vectors of the energy-condition scalars.
    Tf: frame tensor (4,4,N1,N2). Returns dict of 2D arrays (NEC, WEC, SEC, DEC)."""
    dirs = fibonacci_sphere(n_dir)
    eta = np.diag([-1.0, 1, 1, 1])
    T00 = Tf[0, 0]
    T0i = Tf[0, 1:]                     # (3,N1,N2)
    Tij = Tf[1:, 1:]                    # (3,3,N1,N2)
    trT = -T00 + Tij[0, 0] + Tij[1, 1] + Tij[2, 2]
    nec = np.full(T00.shape, np.inf)
    wec = np.full(T00.shape, np.inf)
    sec = np.full(T00.shape, np.inf)
    dec = np.full(T00.shape, np.inf)
    for n in dirs:
        a = np.einsum('i...,i->...', T0i, n)
        b = np.einsum('ij...,i,j->...', Tij, n, n)
        null_val = T00 + 2 * a + b
        nec = np.minimum(nec, null_val)
        # dominant: W^a = -eta^{ac} T_{cb} k^b must be future causal for null k = (1, n)
        k = np.concatenate([[1.0], n])
        Tk = np.einsum('ab...,b->a...', Tf, k)        # T_{ab} k^b
        W = -np.einsum('ac,c...->a...', eta, Tk)      # W^a
        Wnorm = np.einsum('ab,a...,b...->...', eta, W, W)
        W0 = W[0]
        # DEC satisfied where Wnorm <= 0 and W0 >= 0; report the signed margin
        margin = np.where(W0 >= 0, -Wnorm, -np.abs(Wnorm) - np.abs(W0))
        dec = np.minimum(dec, margin)
        for et in etas:
            ch, sh = np.cosh(et), np.sinh(et)
            tval = T00 * ch * ch + 2 * a * ch * sh + b * sh * sh
            wec = np.minimum(wec, tval)
            sec = np.minimum(sec, tval + 0.5 * trT)
    return dict(NEC=nec, WEC=wec, SEC=sec, DEC=dec)


def energy_conditions_eigen(Tf):
    """Hawking-Ellis type I classification: eigen-decompose T^a_b = eta^{ac} T_{cb}.
    Returns rho, p (3), and exact condition margins where type I; NaN elsewhere."""
    eta = np.diag([-1.0, 1, 1, 1])
    N1, N2 = Tf.shape[2:]
    Tmix = np.einsum('ac,cb...->...ab', eta, Tf)       # (N1,N2,4,4)
    w, V = np.linalg.eig(Tmix)
    rho = np.full((N1, N2), np.nan)
    p = np.full((N1, N2, 3), np.nan)
    typeI = np.zeros((N1, N2), bool)
    for i in range(N1):
        for j in range(N2):
            wv = w[i, j]
            if np.max(np.abs(wv.imag)) > 1e-6 * (np.max(np.abs(wv.real)) + 1e-300):
                continue
            wv = wv.real
            Vr = V[i, j].real
            norms = np.einsum('ab,ak,bk->k', eta, Vr, Vr)
            tl = np.where(norms < 0)[0]
            if len(tl) != 1:
                continue
            typeI[i, j] = True
            rho[i, j] = -wv[tl[0]]
            p[i, j] = np.sort(np.delete(wv, tl[0]))
    nec = np.min(rho[..., None] + p, axis=-1)
    wec = np.minimum(rho, nec)
    dec = np.min(rho[..., None] - np.abs(p), axis=-1)
    sec = np.minimum(nec, rho + p.sum(axis=-1))
    return dict(rho=rho, p=p, typeI=typeI, NEC=nec, WEC=wec, SEC=sec, DEC=dec)


# ----------------------------------------------------------------------------
# metric generators (comoving, stationary), on a (Nx,Ny,5) slab
# ----------------------------------------------------------------------------
def grid(xmax, ymax, h, z0=0.0):
    x = np.arange(-xmax, xmax + h / 2, h)
    y = np.arange(-ymax, ymax + h / 2, h)
    z = z0 + h * np.arange(-2, 3)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    return x, y, X, Y, Z


def shape_alcubierre(r, R, sigma):
    return (np.tanh(sigma * (R + r)) + np.tanh(sigma * (R - r))) / (2 * np.tanh(R * sigma))


def minkowski_slab(X):
    g3 = {(0, 0): -np.ones_like(X), (1, 1): np.ones_like(X), (2, 2): np.ones_like(X), (3, 3): np.ones_like(X)}
    return g3


def from_3p1(alpha, beta_lo, gamma):
    """alpha (..), beta_lo list of 3 covariant components, gamma dict (i,j)->array (3x3 symmetric, 1-based 1..3)."""
    shape = alpha.shape
    gam = np.zeros(shape + (3, 3))
    for (i, j), v in gamma.items():
        gam[..., i - 1, j - 1] = v
        gam[..., j - 1, i - 1] = v
    gaminv = np.linalg.inv(gam)
    bl = np.stack(beta_lo, axis=-1)
    bu = np.einsum('...ij,...j->...i', gaminv, bl)
    g3 = {(0, 0): -alpha ** 2 + np.einsum('...i,...i->...', bu, bl)}
    for i in range(3):
        g3[(0, i + 1)] = bl[..., i]
    for (i, j), v in gamma.items():
        g3[(i, j)] = v
    return g3


def alcubierre_comoving(X, Y, Z, v, R, sigma):
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    f = shape_alcubierre(r, R, sigma)
    beta_x = v * (1 - f)
    one = np.ones_like(X)
    return from_3p1(one, [beta_x, 0 * X, 0 * X], {(1, 1): one, (2, 2): one, (3, 3): one})


def vandenbroeck_comoving(X, Y, Z, v, R1, sigma1, R2, sigma2, A, as_coded=False):
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    B = 1 + A * shape_alcubierre(r, R1, sigma1)
    f = shape_alcubierre(r, R2, sigma2)
    g3 = {(1, 1): B ** 2, (2, 2): B ** 2, (3, 3): B ** 2, (0, 1): B ** 2 * v * (1 - f)}
    if as_coded:
        g3[(0, 0)] = -(1 - B ** 2 * (v * f) ** 2)          # repository main branch, Oct 2026
    else:
        g3[(0, 0)] = -1 + B ** 2 * v ** 2 * (1 - f) ** 2    # Galilean comoving form of Van Den Broeck's line element
    return g3


def modified_time_comoving(X, Y, Z, v, R, sigma, A, as_coded=False):
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    f = shape_alcubierre(r, R, sigma)
    Nlapse = (1 - f) + f / A
    one = np.ones_like(X)
    g3 = {(1, 1): one, (2, 2): one, (3, 3): one, (0, 1): v * (1 - f)}
    if as_coded:
        g3[(0, 0)] = -Nlapse ** 2 + (f * v) ** 2
    else:
        g3[(0, 0)] = -Nlapse ** 2 + (v * (1 - f)) ** 2
    return g3


def lentz_template(x, y, s):
    """Vectorised port of getWarpFactorByRegion (Warp Factory, Metrics/Lentz)."""
    ya = np.abs(y)
    WFX = np.zeros_like(x)
    WFY = np.zeros_like(x)
    done = np.zeros_like(x, bool)
    conds = [
        ((x >= s) & (x <= 2 * s) & (x - s >= ya), -2.0, 0.0),
        ((x > s) & (x <= 2 * s) & (x - s <= ya) & (-ya + 3 * s >= x), -1.0, 1.0),
        ((x > 0) & (x <= s) & (x + s > ya) & (-ya + s < x), 0.0, 1.0),
        ((x > 0) & (x <= s) & (x + s <= ya) & (-ya + 3 * s >= x), -0.5, 0.5),
        ((x > -s) & (x <= 0) & (-x + s < ya) & (-ya + 3 * s >= -x), 0.5, 0.5),
        ((x > -s) & (x <= 0) & (x + s <= ya) & (-ya + s >= x), 1.0, 0.0),
        ((x >= -s) & (x <= s) & (x + s > ya), 1.0, 0.0),
    ]
    for cnd, wx, wy in conds:
        m = cnd & ~done
        WFX[m] = wx
        WFY[m] = wy
        done |= m
    return WFX, np.sign(y) * WFY


def lentz_comoving(X, Y, Z, v, s, smooth_cells=0, h=None):
    WFX, WFY = lentz_template(X, Y, s)
    beta_x = v * (1 - WFX)
    beta_y = v * WFY
    if smooth_cells > 0:
        from scipy.ndimage import gaussian_filter
        beta_x = gaussian_filter(beta_x, sigma=(smooth_cells, smooth_cells, 0), mode='nearest')
        beta_y = gaussian_filter(beta_y, sigma=(smooth_cells, smooth_cells, 0), mode='nearest')
    one = np.ones_like(X)
    return from_3p1(one, [beta_x, beta_y, 0 * X], {(1, 1): one, (2, 2): one, (3, 3): one})


def matlab_smooth(y, span):
    """MATLAB smooth(y,span,'moving'): centred moving average, odd span, shrinking symmetric windows at the ends."""
    span = int(np.floor(span))
    if span % 2 == 0:
        span -= 1
    span = max(span, 1)
    n = len(y)
    cs = np.concatenate([[0.0], np.cumsum(y)])
    out = np.empty(n)
    hw = span // 2
    idx = np.arange(n)
    w = np.minimum(hw, np.minimum(idx, n - 1 - idx))
    lo = idx - w
    hi = idx + w
    out = (cs[hi + 1] - cs[lo]) / (hi - lo + 1)
    return out


def compact_sigmoid(r, R1, R2, sigma, Rbuff):
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        inner = (r > R1 + Rbuff) & (r < R2 - Rbuff)
        arg = ((R2 - R1 - 2 * Rbuff) * (sigma + 2)) / 2 * (1.0 / (r - R2 + Rbuff) + 1.0 / (r - R1 - Rbuff))
        sg = 1.0 / (np.exp(arg) + 1.0)
        f = np.abs(np.where(inner, sg, 0.0) + (r >= R2 - Rbuff) - 1.0)
    return f


def warpshell_profiles(m, R1, R2, Rbuff, sigma, smoothFactor, rmax, n=10 ** 5):
    """Port of the radial construction in metricGet_WarpShellComoving."""
    rs = np.linspace(0, rmax, n)
    rho = (m / (4 / 3 * np.pi * (R2 ** 3 - R1 ** 3))) * ((rs > R1) & (rs < R2))
    from scipy.integrate import cumulative_trapezoid
    M = cumulative_trapezoid(4 * np.pi * rho * rs ** 2, rs, initial=0)
    Mend = M[-1]
    Rr = R2
    with np.errstate(invalid='ignore'):
        P = C ** 2 * rho * ((Rr * np.sqrt(Rr - 2 * GN * Mend / C ** 2) - np.sqrt(Rr ** 3 - 2 * GN * Mend * rs ** 2 / C ** 2))
                            / (np.sqrt(Rr ** 3 - 2 * GN * Mend * rs ** 2 / C ** 2) - 3 * Rr * np.sqrt(Rr - 2 * GN * Mend / C ** 2))) * (rs < Rr)
    P = np.nan_to_num(P)
    for _ in range(4):
        rho = matlab_smooth(rho, 1.79 * smoothFactor)
    for _ in range(4):
        P = matlab_smooth(P, smoothFactor)
    M = cumulative_trapezoid(4 * np.pi * rho * rs ** 2, rs, initial=0)
    M[M < 0] = M.max()
    S = compact_sigmoid(rs, R1, R2, sigma, Rbuff)
    S = matlab_smooth(matlab_smooth(S, smoothFactor), smoothFactor)
    with np.errstate(divide='ignore', invalid='ignore'):
        B = 1.0 / (1 - 2 * GN * M / rs / C ** 2)
    B[0] = 1.0
    with np.errstate(divide='ignore', invalid='ignore'):
        dalpha = (GN * M / C ** 2 + 4 * np.pi * GN * rs ** 3 * P / C ** 4) / (rs ** 2 - 2 * GN * M * rs / C ** 2)
    dalpha[0] = 0.0
    aT = cumulative_trapezoid(dalpha, rs, initial=0)
    Cc = 0.5 * np.log(1 - 2 * GN * Mend / rs[-1] / C ** 2)
    a = aT + (Cc - aT[-1])
    A = -np.exp(2 * a)
    return dict(r=rs, rho=rho, P=P, M=M, S=S, A=A, B=B)


def warpshell_comoving(X, Y, Z, prof, vWarp, doWarp=True):
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    A = np.interp(r, prof['r'], prof['A'])
    B = np.interp(r, prof['r'], prof['B'])
    S = np.interp(r, prof['r'], prof['S'])
    rr = np.where(r > 0, r, 1.0)
    nx, ny, nz = X / rr, Y / rr, Z / rr
    g3 = {(0, 0): A,
          (1, 1): 1 + (B - 1) * nx * nx, (2, 2): 1 + (B - 1) * ny * ny, (3, 3): 1 + (B - 1) * nz * nz,
          (1, 2): (B - 1) * nx * ny, (1, 3): (B - 1) * nx * nz, (2, 3): (B - 1) * ny * nz}
    if doWarp:
        g3[(0, 1)] = -S * vWarp
    return g3


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------
def analyse(g3, hx, hy, hz, margin=4, n_dir=400, eigen=True):
    g, dg, d2g = slab_derivs(g3, hx, hy, hz)
    ginv, Ric, Rs, Gmn = einstein_from_derivs(g, dg, d2g)
    T = KAPPA * Gmn
    e, alpha, beta_up = eulerian_frame(g)
    Tf = frame_tensor(T, e)
    sl = (slice(margin, -margin), slice(margin, -margin))
    Tf_c = Tf[(slice(None), slice(None)) + sl]
    out = dict(g=g, T=T, Tf=Tf, alpha=alpha, beta_up=beta_up, slice=sl)
    out['rho_E'] = Tf[0, 0]
    out['sampled'] = energy_conditions_sampled(Tf_c, n_dir=n_dir)
    if eigen:
        out['eigen'] = energy_conditions_eigen(Tf_c)
    return out


def summary(out, label):
    sl = out['slice']
    rho = out['rho_E'][sl]
    s = out['sampled']
    line = {
        'metric': label,
        'rho_E min [J/m^3]': float(np.nanmin(rho)),
        'rho_E max [J/m^3]': float(np.nanmax(rho)),
        'NEC min': float(np.nanmin(s['NEC'])),
        'WEC min': float(np.nanmin(s['WEC'])),
        'SEC min': float(np.nanmin(s['SEC'])),
        'DEC min': float(np.nanmin(s['DEC'])),
    }
    if 'eigen' in out:
        ev = out['eigen']
        line['typeI fraction'] = float(np.mean(ev['typeI']))
        for k in ('NEC', 'WEC', 'SEC', 'DEC'):
            line[k + ' min (eigen)'] = float(np.nanmin(ev[k]))
    return line
