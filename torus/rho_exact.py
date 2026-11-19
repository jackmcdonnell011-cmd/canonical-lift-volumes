"""Exact (non-numerical) Ghys map rho : pi_1(S^3 - (trefoil u K_w)) -> PSL2(Z).

rho kills the meridian of K, so it factors through the filled group pi_1(S^3 - T).  SnapPy fills the
K cusp on the same triangulation, simplifies that group (normally to <a, b | a^2 b^3>) and records how
every original (face-dual) generator is written in the simplified ones.  The unfilled group's
generators are also known as words in the original generators.  So:
  1. send the filled generators to PSL2(Z): an element of order 2 and one of order 3 that generate
     PSL2(Z) = Z/2 * Z/3, i.e. (S, R) or (S, R^2) up to conjugation, R = S X;
  2. conjugate so that the meridian of T goes to X (use the automorphism J if it goes to X^-1);
  3. carry this back to the unfilled generators through the original generators.
The result is checked exactly by rho.verify (every relator -> +-I, mu_K -> +-I, mu_T -> +-X), and
ModularLink then checks that lambda_K is the word, so a wrong answer cannot get through.
"""
from math import gcd
import psl2

R = psl2.mul(psl2.S, psl2.X)


def _word_int(word, mats):
    M = psl2.I2
    for i in word:
        A = mats[abs(i) - 1]
        M = psl2.mul(M, A if i > 0 else psl2.inv(A))
    return M


def _to_infinity(P):
    """h in PSL2(Z) with h P h^-1 = +-[[1, k], [0, 1]], for a parabolic P."""
    a, b, c, d = P
    if a + d < 0:
        a, b, c, d = -a, -b, -c, -d
    if c == 0:
        return psl2.I2
    # fixed point (a - d) / (2c) = p/q in lowest terms
    p, q = a - d, 2 * c
    g = gcd(p, q)
    p, q = p // g, q // g
    # s p - r q = 1
    def egcd(x, y):
        if y == 0:
            return (x, 1, 0)
        g, s, t = egcd(y, x % y)
        return (g, t, s - (x // y) * t)
    g, s, t = egcd(p, q)           # s p + t q = g = +-1
    s, t = s * g, t * g
    h = (s, t, -q, p)               # det = s p + t q = 1, h(p/q) = oo
    return h


def find_rho_exact(M, cuspK=0, cuspT=1):
    """M: the link exterior (unfilled, cusps K and T).  Returns gm {generator: matrix} for
    M.fundamental_group(), or None."""
    G = M.fundamental_group()
    Mf = M.copy()
    Mf.dehn_fill((1, 0), cuspK)
    F = Mf.fundamental_group()
    fg = F.generators()
    orig_in_F = F.original_generators(as_int_list=True)
    gens_in_orig = G.generators_in_originals(as_int_list=True)
    muT = G.peripheral_curves()[cuspT][0]
    if len(fg) != 2:
        return None
    import rho as rhomod
    for order2 in (0, 1):
        for e in (1, 2):
            fmats = [None, None]
            fmats[order2] = psl2.S
            fmats[1 - order2] = R if e == 1 else psl2.mul(R, R)
            fm = dict(zip(fg, fmats))
            if not all(psl2.eq_pm(psl2.word_matrix(r, fm), psl2.I2) for r in F.relators()):
                continue
            orig = [_word_int(wd, fmats) for wd in orig_in_F]
            gm = {g: _word_int(wd, orig) for g, wd in zip(G.generators(), gens_in_orig)}
            # normalise with the unfilled group's own meridian of T (the filled group's is another conjugate)
            P = psl2.word_matrix(muT, gm)
            if abs(psl2.tr(P)) != 2 or psl2.eq_pm(P, psl2.I2):
                continue
            h = _to_infinity(P)
            hi = psl2.inv(h)
            Q = psl2.mul(psl2.mul(h, P), hi)
            k = Q[1] if Q[0] > 0 else -Q[1]
            if abs(k) != 1:
                continue
            gm = {g: psl2.mul(psl2.mul(h, A), hi) for g, A in gm.items()}
            if k == -1:
                gm = {g: psl2.conjJ(A) for g, A in gm.items()}
            if rhomod.verify(G, gm, cuspK, cuspT):
                return gm
    return None
