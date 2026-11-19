"""Find the 'Ghys' homomorphism rho : pi_1(S^3 - (trefoil u K_w)) -> PSL2(Z).

Theory used: pi_1(S^3 - T) = B_3, and B_3 / centre = PSL2(Z).  Killing the
meridian of K_w gives pi_1(S^3 - T).  Up to conjugation in PGL2(C) there is a
unique irreducible SL2(C) rep of B_3 (mod sign) in which a meridian of T is a
nontrivial parabolic, and it is the standard one onto PSL2(Z).  So we solve
numerically for a rep of SnapPy's presentation with
    rho(mu_K) = +-I,  rho(mu_T) = +-X,  relators -> +-I,  rep irreducible,
then conjugate by a unipotent to make it integral and verify EXACTLY.
Any isomorphism G/<<mu_K>>/Z -> PSL2(Z) differs from the true Ghys one by an
automorphism of PSL2(Z) (inner, or inner o J); see notes in README section of
the report for why that does not affect covers once we label by rho(lambda_K).
"""
import itertools
import numpy as np
from scipy.optimize import least_squares
from fractions import Fraction
import psl2


def _wmat(word, mats):
    M = np.eye(2)
    for ch in word:
        i = ord(ch.lower()) - 97
        A = mats[i]
        if ch.isupper():
            A = np.array([[A[1, 1], -A[0, 1]], [-A[1, 0], A[0, 0]]])
        M = M @ A
    return M


def _compile(word):
    return [(ord(ch.lower()) - 97, ch.isupper()) for ch in word]


def _pw(cw, mats):
    a, b, c, d = 1.0, 0.0, 0.0, 1.0
    for i, invflag in cw:
        e, f, g, h = mats[i]
        if invflag:
            e, f, g, h = h, -f, -g, e
        a, b, c, d = a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h
    return a, b, c, d


def find_rho(G, cuspK=0, cuspT=1, tries=2000, seed=0, stats=None, max_nfev=150):
    gens = G.generators()
    k = len(gens)
    rels = G.relators()
    per = G.peripheral_curves()
    muK = per[cuspK][0]
    muT = per[cuspT][0]
    rng = np.random.default_rng(seed)
    Xn = np.array([[1., 1.], [0., 1.]])
    E = np.eye(2)

    crels = [_compile(r) for r in rels]
    cK, cT = _compile(muK), _compile(muT)

    def resid(v):
        mats = [tuple(v[4 * i:4 * i + 4]) for i in range(k)]
        out = [m[0] * m[3] - m[1] * m[2] - 1 for m in mats]
        for cw in crels + [cK]:
            a, b, c, d = _pw(cw, mats)
            s_ = 1.0 if a + d >= 0 else -1.0
            out += [a - s_, b, c, d - s_]
        a, b, c, d = _pw(cT, mats)
        s_ = 1.0 if a + d >= 0 else -1.0
        out += [a - s_, b - s_, c, d - s_]
        return np.array(out)

    for t in range(tries):
        scale = 3.0
        v0 = rng.uniform(-scale, scale, 4 * k)
        try:
            sol = least_squares(resid, v0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=max_nfev)
        except Exception:
            continue
        if np.max(np.abs(sol.fun)) > 1e-9:
            continue
        mats = [sol.x[4 * i:4 * i + 4].reshape(2, 2) for i in range(k)]
        irr = False
        for i in range(k):
            for j in range(k):
                C = mats[i] @ mats[j] @ np.linalg.inv(mats[i]) @ np.linalg.inv(mats[j])
                if abs(np.trace(C) - 2) > 1e-6:
                    irr = True
        if not irr:
            continue
        ints = _integerize(mats)
        if ints is None:
            continue
        gm = {chr(97 + i): ints[i] for i in range(k)}
        if verify(G, gm, cuspK, cuspT):
            if stats is not None:
                stats['tries'] = t + 1
            return gm
    return None


def _integerize(mats):
    """mats are U M_i U^{-1} with U=[[1,s],[0,1]], M_i integral. Recover M_i."""
    cands = []
    for m in mats:
        c = m[1, 0]
        if abs(c) > 0.5:
            cands.append(m)
    if not cands:
        return None
    m = min(cands, key=lambda m: abs(m[1, 0]))
    c = int(round(m[1, 0]))
    a1 = m[0, 0]
    # a1 = a + s c  ->  s = (a1 - a)/c for integer a
    base = int(np.floor(a1))
    for a in range(base - abs(c) - 1, base + abs(c) + 2):
        s = (a1 - a) / c
        U = np.array([[1, s], [0, 1]])
        Ui = np.array([[1, -s], [0, 1]])
        out = []
        ok = True
        for mm in mats:
            r = Ui @ mm @ U
            ri = np.round(r)
            if np.max(np.abs(r - ri)) > 1e-6:
                ok = False
                break
            out.append(tuple(int(x) for x in ri.ravel()))
        if ok:
            return out
    return None


def verify(G, gm, cuspK=0, cuspT=1):
    for g, A in gm.items():
        if A[0] * A[3] - A[1] * A[2] != 1:
            return False
    for r in G.relators():
        if not psl2.eq_pm(psl2.word_matrix(r, gm), psl2.I2):
            return False
    per = G.peripheral_curves()
    if not psl2.eq_pm(psl2.word_matrix(per[cuspK][0], gm), psl2.I2):
        return False
    if not psl2.eq_pm(psl2.word_matrix(per[cuspT][0], gm), psl2.X):
        return False
    return True


def twist(gm):
    """Compose with the outer automorphism J of PSL2(Z)."""
    return {g: psl2.conjJ(A) for g, A in gm.items()}
