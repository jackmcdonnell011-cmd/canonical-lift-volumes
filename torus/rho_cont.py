"""Fallback solver for rho: path-following along the character curve of the trefoil group.

Irreducible SL2(C) reps of B_3 = pi_1(S^3 - T) form one curve parametrised by
t = tr(mu_T).  We solve at t = 2cos(theta0) (elliptic mu_T, bounded solutions, easy to
hit from random starts), then follow t -> 2 (theta -> 0) by Newton steps, normalising
    rho(mu_T) = [[e^{i theta}, 1], [0, e^{-i theta}]]   -> X at theta = 0.
At theta = 0 the result is a conjugate of the integral rep; integerise as before.
"""
import numpy as np
from scipy.optimize import least_squares
from rho import _compile, verify
import cmath
import math


def _pw(cw, mats):
    a, b, c, d = 1, 0, 0, 1
    for i, inv in cw:
        e, f, g, h = mats[i]
        if inv:
            e, f, g, h = h, -f, -g, e
        a, b, c, d = a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h
    return a, b, c, d


def _make(G, cuspK, cuspT):
    k = len(G.generators())
    per = G.peripheral_curves()
    crels = [_compile(r) for r in G.relators()]
    cK, cT = _compile(per[cuspK][0]), _compile(per[cuspT][0])

    def resid(v, theta):
        # theta may be complex: lam = exp(i*theta)
        z = v[:4 * k] + 1j * v[4 * k:]
        mats = [tuple(z[4 * i:4 * i + 4]) for i in range(k)]
        out = [m[0] * m[3] - m[1] * m[2] - 1 for m in mats]
        for cw in crels + [cK]:
            a, b, c, d = _pw(cw, mats)
            s_ = 1 if (a + d).real >= 0 else -1
            out += [a - s_, b, c, d - s_]
        a, b, c, d = _pw(cT, mats)
        lam = cmath.exp(1j * theta)
        s_ = 1 if ((a + d) * complex(lam + 1 / lam).conjugate()).real >= 0 else -1
        out += [a - s_ * lam, b - s_, c, d - s_ / lam]
        out = np.array(out)
        return np.concatenate([out.real, out.imag])
    return k, resid


def _irreducible(mats):
    k = len(mats)
    for i in range(k):
        for j in range(k):
            A, B = mats[i], mats[j]
            C = A @ B @ np.linalg.inv(A) @ np.linalg.inv(B)
            if abs(np.trace(C) - 2) > 1e-6:
                return True
    return False


def find_rho_cont(G, cuspK=0, cuspT=1, theta0=np.pi / 3, tries=200, steps=80, seed=0, detour=0.3, stats=None):
    rng = np.random.default_rng(seed)
    k, resid = _make(G, cuspK, cuspT)
    for t in range(tries):
        v0 = rng.uniform(-1.5, 1.5, 8 * k)
        sol = least_squares(lambda v: resid(v, theta0), v0, method='lm', max_nfev=200,
                            xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if np.max(np.abs(sol.fun)) > 1e-9:
            continue
        z = sol.x[:4 * k] + 1j * sol.x[4 * k:]
        if not _irreducible([z[4 * i:4 * i + 4].reshape(2, 2) for i in range(k)]):
            continue
        v = sol.x
        ok = True
        # adaptive path-following in the complex theta-plane: detour away from the
        # reducible/irreducible branch point (theta = pi/6) of the trefoil character curve
        path = lambda sp: theta0 * (1 - sp) - 1j * detour * math.sin(math.pi * sp)
        sp, h, nsteps = 0.0, 1.0 / steps, 0
        while sp < 1 - 1e-12 and ok:
            h = min(h, 1 - sp)
            s2 = least_squares(lambda u: resid(u, path(sp + h)), v, method='lm', max_nfev=200,
                               xtol=1e-15, ftol=1e-15, gtol=1e-15)
            nsteps += 1
            if np.max(np.abs(s2.fun)) < 1e-9 and np.max(np.abs(s2.x - v)) < 2.0 + 0.5 * np.max(np.abs(v)):
                v, sp = s2.x, sp + h
                h = min(h * 1.5, 0.05)
            else:
                h /= 2
                if h < 1e-4 or nsteps > 400:
                    ok = False
        if not ok:
            continue
        z = v[:4 * k] + 1j * v[4 * k:]
        mats = [z[4 * i:4 * i + 4].reshape(2, 2) for i in range(k)]
        ints = _integerize_c(mats)
        if ints is None:
            continue
        gm = {chr(97 + i): ints[i] for i in range(k)}
        if verify(G, gm, cuspK, cuspT):
            if stats is not None:
                stats['tries'] = t + 1
            return gm
    return None


def _integerize_c(mats):
    cands = [m for m in mats if abs(m[1, 0]) > 0.5]
    if not cands:
        return None
    m = min(cands, key=lambda m: abs(m[1, 0]))
    c = m[1, 0]
    if abs(c - round(c.real)) > 1e-5:
        return None
    c = int(round(c.real))
    a1 = m[0, 0]
    base = int(np.floor(a1.real))
    for a in range(base - abs(c) - 1, base + abs(c) + 2):
        s = (a1 - a) / c
        U = np.array([[1, s], [0, 1]])
        Ui = np.array([[1, -s], [0, 1]])
        out = []
        for mm in mats:
            r = Ui @ mm @ U
            ri = np.round(r.real)
            if np.max(np.abs(r - ri)) > 1e-6:
                out = None
                break
            out.append(tuple(int(x) for x in ri.ravel()))
        if out:
            return out
    return None
