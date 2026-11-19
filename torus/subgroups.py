"""Finite-index subgroups Gamma of PSL2(Z), given by the right action of PSL2(Z)
on the cosets Gamma\\PSL2(Z).  Each subgroup object exposes
    .d                   index
    .act(i, M)           coset index i acted on (on the right) by matrix M
    .perm(M)             permutation list p with p[i] = i.M
    .signature()         (genus, cusps, e2, e3) of H/Gamma
"""
from math import gcd
import psl2


class CosetAction:
    def __init__(self, name, base, right_act):
        """base: hashable point; right_act(pt, M) -> pt (must be a right action of PSL2(Z))."""
        self.name = name
        self._ra = right_act
        pts = [base]
        reps = [psl2.I2]   # reps[i] = h with (coset 0).h = coset i, i.e. coset i = Gamma h
        idx = {base: 0}
        k = 0
        while k < len(pts):
            p = pts[k]
            for g in (psl2.S, psl2.X, psl2.inv(psl2.X)):
                q = right_act(p, g)
                if q not in idx:
                    idx[q] = len(pts)
                    pts.append(q)
                    reps.append(psl2.mul(reps[k], g))
            k += 1
        self.pts, self.idx, self.d, self.reps = pts, idx, len(pts), reps
        self.pS = self.perm(psl2.S)
        self.pT = self.perm(psl2.X)

    def act(self, i, M):
        return self.idx[self._ra(self.pts[i], M)]

    def perm(self, M):
        return [self.act(i, M) for i in range(len(self.pts))]

    def signature(self):
        d = self.d
        pS, pT = self.pS, self.pT
        pU = [pT[pS[i]] for i in range(d)]  # i.(S T) : first S then T  (U = S T has order 3)
        e2 = sum(1 for i in range(d) if pS[i] == i)
        e3 = sum(1 for i in range(d) if pU[i] == i)
        # cusps = cycles of T
        seen, cusps, widths = set(), 0, []
        for i in range(d):
            if i in seen:
                continue
            j, w = i, 0
            while j not in seen:
                seen.add(j)
                j = pT[j]
                w += 1
            cusps += 1
            widths.append(w)
        # Riemann-Hurwitz: g = 1 + d/12 - e2/4 - e3/3 - cusps/2
        g12 = 12 + d - 3 * e2 - 4 * e3 - 6 * cusps
        assert g12 % 12 == 0, (d, e2, e3, cusps)
        return dict(genus=g12 // 12, cusps=cusps, e2=e2, e3=e3, widths=sorted(widths), index=d)


def _proj_mod(c, d, N):
    """Normalise (c:d) in P^1(Z/N) by choosing minimal representative under units."""
    best = None
    for u in range(1, N):
        if gcd(u, N) == 1:
            r = ((u * c) % N, (u * d) % N)
            if best is None or r < best:
                best = r
    return best


def Gamma0(N):
    def ra(p, M):
        c, d = p
        a, b, cc, dd = M
        return _proj_mod(c * a + d * cc, c * b + d * dd, N)
    return CosetAction(f'Gamma0({N})', _proj_mod(0, 1, N), ra)


def Gamma1(N):
    def norm(c, d):
        r1 = (c % N, d % N)
        r2 = ((-c) % N, (-d) % N)
        return min(r1, r2)

    def ra(p, M):
        c, d = p
        a, b, cc, dd = M
        return norm(c * a + d * cc, c * b + d * dd)
    return CosetAction(f'Gamma1({N})', norm(0, 1), ra)


def GammaN(N):
    def norm(A):
        r1 = tuple(x % N for x in A)
        r2 = tuple((-x) % N for x in A)
        return min(r1, r2)

    def ra(p, M):
        return norm(psl2.mul(p, M))
    return CosetAction(f'Gamma({N})', norm(psl2.I2), ra)


def chi6(M):
    """The abelianisation PSL2(Z) -> Z/6 with S -> 3, T -> 5 (U=ST -> 2)."""
    s = 0
    for t, k in psl2.st_word(M):
        s += 3 if t == 'S' else 5 * k
    return s % 6


def Commutator():
    return CosetAction("Gamma'", 0, lambda p, M: (p + chi6(M)) % 6)


def Intersect(A, B):
    """Gamma_A cap Gamma_B: action on orbit of (base,base) in the product."""
    def ra(p, M):
        return (A._ra(p[0], M), B._ra(p[1], M))
    return CosetAction(f'{A.name} cap {B.name}', (A.pts[0], B.pts[0]), ra)


def FromPerms(name, pS, pT):
    """Arbitrary subgroup from permutations of S and T (right action, p[i] = i.g)."""
    pS, pT = list(pS), list(pT)
    d = len(pS)
    pTi = [0] * d
    for i, j in enumerate(pT):
        pTi[j] = i

    def ra(p, M):
        for t, k in psl2.st_word(M):
            if t == 'S':
                p = pS[p]
            else:
                pp = pT if k > 0 else pTi
                for _ in range(abs(k)):
                    p = pp[p]
        return p
    return CosetAction(name, 0, ra)
