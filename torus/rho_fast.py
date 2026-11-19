"""Fast, robust solver for the Ghys map rho : pi_1(S^3 - (T u K_w)) -> PSL2(Z).

Same mathematics as rho.py / rho_cont.py (see there), but with
  * an analytic (holomorphic) Jacobian of all word maps, and a damped complex
    Gauss-Newton iteration instead of scipy's finite-difference LM;
  * path-following along the trefoil character curve, parametrised by
    rho(mu_T) = [[lam, 1], [0, 1/lam]], lam = exp(i*theta), from an elliptic start
    theta0 = pi/3 to theta = 0, detouring into complex theta round the branch point
    theta = pi/6 where the irreducible and reducible characters meet;
  * a cheap direct attempt at theta = 0 first.
The result is rounded to integers and verified EXACTLY (rho.verify), so numerical
error can only make the solver fail, never return a wrong map.
"""
import math
import cmath
import numpy as np
from rho import _compile, verify


def _mat_ops():
    # derivative of inv(A) = [[d,-b],[-c,a]] w.r.t. (a,b,c,d) (valid on det=1)
    E = [np.array([[1, 0], [0, 0]], complex), np.array([[0, 1], [0, 0]], complex),
         np.array([[0, 0], [1, 0]], complex), np.array([[0, 0], [0, 1]], complex)]
    Einv = [np.array([[0, 0], [0, 1]], complex), np.array([[0, -1], [0, 0]], complex),
            np.array([[0, 0], [-1, 0]], complex), np.array([[1, 0], [0, 0]], complex)]
    return E, Einv


E_, EINV_ = _mat_ops()


def _word_and_jac(cw, A, Ainv, k):
    """Value (2x2) and derivative (4k x 2 x 2) of the word map."""
    L = len(cw)
    pre = [np.eye(2, dtype=complex)]
    for (i, inv) in cw:
        pre.append(pre[-1] @ (Ainv[i] if inv else A[i]))
    suf = [np.eye(2, dtype=complex)] * (L + 1)
    for t in range(L - 1, -1, -1):
        i, inv = cw[t]
        suf[t] = (Ainv[i] if inv else A[i]) @ suf[t + 1]
    J = np.zeros((4 * k, 2, 2), complex)
    for t, (i, inv) in enumerate(cw):
        P, Sx = pre[t], suf[t + 1]
        Es = EINV_ if inv else E_
        for e in range(4):
            J[4 * i + e] += P @ Es[e] @ Sx
    return pre[-1], J


class System:
    def __init__(self, G, cuspK=0, cuspT=1):
        self.G = G
        self.k = len(G.generators())
        per = G.peripheral_curves()
        self.rels = [_compile(r) for r in G.relators()] + [_compile(per[cuspK][0])]
        self.cT = _compile(per[cuspT][0])
        self.cuspK, self.cuspT = cuspK, cuspT

    def F(self, z, theta):
        k = self.k
        A = [z[4 * i:4 * i + 4].reshape(2, 2) for i in range(k)]
        Ainv = [np.array([[a[1, 1], -a[0, 1]], [-a[1, 0], a[0, 0]]]) for a in A]
        res, jac = [], []
        for i in range(k):  # det = 1
            a, b, c, d = z[4 * i:4 * i + 4]
            res.append(a * d - b * c - 1)
            row = np.zeros(4 * k, complex)
            row[4 * i:4 * i + 4] = [d, -c, -b, a]
            jac.append(row)
        for cw in self.rels:
            W, J = _word_and_jac(cw, A, Ainv, k)
            s = 1 if (W[0, 0] + W[1, 1]).real >= 0 else -1
            T = W - s * np.eye(2)
            res.extend(T.ravel())
            jac.extend(J.reshape(4 * k, 4).T)
        lam = cmath.exp(1j * theta)
        W, J = _word_and_jac(self.cT, A, Ainv, k)
        target = np.array([[lam, 1], [0, 1 / lam]])
        s = 1 if ((W[0, 0] + W[1, 1]) * np.conj(lam + 1 / lam)).real >= 0 else -1
        res.extend((W - s * target).ravel())
        jac.extend(J.reshape(4 * k, 4).T)
        return np.array(res), np.array(jac)

    def newton(self, z, theta, maxit=40, tol=1e-11):
        r, J = self.F(z, theta)
        nr = np.linalg.norm(r)
        for _ in range(maxit):
            if nr < tol:
                return z, True
            dz = np.linalg.lstsq(J, -r, rcond=None)[0]
            t = 1.0
            while t > 1e-4:
                z2 = z + t * dz
                r2, J2 = self.F(z2, theta)
                n2 = np.linalg.norm(r2)
                if n2 < nr:
                    z, r, J, nr = z2, r2, J2, n2
                    break
                t /= 2
            else:
                return z, False
        return z, nr < tol

    def irreducible(self, z):
        k = self.k
        A = [z[4 * i:4 * i + 4].reshape(2, 2) for i in range(k)]
        for a in A:
            for b in A:
                C = a @ b @ np.linalg.inv(a) @ np.linalg.inv(b)
                if abs(np.trace(C) - 2) > 1e-6:
                    return True
        return False


def _integerize(z, k):
    mats = [z[4 * i:4 * i + 4].reshape(2, 2) for i in range(k)]
    cands = [m for m in mats if abs(m[1, 0]) > 0.5]
    if not cands:
        return None
    m = min(cands, key=lambda m: abs(m[1, 0]))
    c = m[1, 0]
    if abs(c - round(c.real)) > 1e-6:
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


def find_rho_fast(G, cuspK=0, cuspT=1, seed=0, direct_tries=40, cont_tries=60,
                  theta0=math.pi / 3, detour=0.3, stats=None, max_steps=600):
    rng = np.random.default_rng(seed)
    sysm = System(G, cuspK, cuspT)
    k = sysm.k

    def finish(z, how, t):
        ints = _integerize(z, k)
        if ints is None:
            return None
        gm = {chr(97 + i): ints[i] for i in range(k)}
        if verify(G, gm, cuspK, cuspT):
            if stats is not None:
                stats.update(method=how, tries=t + 1)
            return gm
        return None

    # 1. direct attempts at theta = 0 from real random starts
    for t in range(direct_tries):
        z0 = rng.uniform(-3, 3, 4 * k).astype(complex)
        z, ok = sysm.newton(z0, 0.0, maxit=60)
        if ok and sysm.irreducible(z):
            gm = finish(z, 'direct', t)
            if gm:
                return gm
    # 2. continuation from an elliptic start
    path = lambda sp: theta0 * (1 - sp) - 1j * detour * math.sin(math.pi * sp)
    for t in range(cont_tries):
        z0 = rng.uniform(-1.5, 1.5, 4 * k) + 1j * rng.uniform(-1.5, 1.5, 4 * k)
        z, ok = sysm.newton(z0, path(0.0), maxit=60)
        if not (ok and sysm.irreducible(z)):
            continue
        sp, h, steps, zprev, hprev = 0.0, 0.02, 0, None, None
        alive = True
        while sp < 1 - 1e-12 and alive:
            h = min(h, 1 - sp)
            pred = z if zprev is None else z + (z - zprev) * (h / hprev)
            z2, ok = sysm.newton(pred, path(sp + h), maxit=12)
            steps += 1
            if ok and np.max(np.abs(z2 - z)) < 1.0 + 0.3 * np.max(np.abs(z)):
                zprev, hprev = z, h
                z, sp = z2, sp + h
                h = min(h * 1.6, 0.08)
            else:
                h /= 2
                zprev = None
                if h < 1e-5 or steps > max_steps:
                    alive = False
        if alive:
            gm = finish(z, 'continuation', t)
            if gm:
                return gm
    return None
