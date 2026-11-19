"""Exact PSL2(Z) helpers: matrices as 4-tuples (a,b,c,d) of Python ints."""
from math import isqrt

I2 = (1, 0, 0, 1)
X = (1, 1, 0, 1)
Y = (1, 0, 1, 1)
S = (0, -1, 1, 0)
J = (1, 0, 0, -1)  # in PGL2(Z): z -> -z (orientation reversing)


def mul(A, B):
    a, b, c, d = A
    e, f, g, h = B
    return (a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h)


def inv(A):
    a, b, c, d = A
    return (d, -b, -c, a)


def neg(A):
    return tuple(-x for x in A)


def tr(A):
    return A[0] + A[3]


def eq_pm(A, B):
    return A == B or A == neg(B)


def word_matrix(w, gens):
    """w: string over letters; lowercase = generator, uppercase = inverse."""
    M = I2
    for ch in w:
        if ch.islower():
            M = mul(M, gens[ch])
        else:
            M = mul(M, inv(gens[ch.lower()]))
    return M


def xy_matrix(w):
    M = I2
    for ch in w:
        M = mul(M, X if ch == 'X' else Y)
    return M


def conjJ(A):
    a, b, c, d = A
    return (a, -b, -c, d)


# ---------- decomposition into S and T=X (for permutation actions) ----------
def st_word(M):
    """Return list of tokens ('T',k) / ('S',1) with M = product (in PSL2(Z))."""
    a, b, c, d = M
    toks = []
    while c != 0:
        q = a // c
        # left-multiply by T^{-q}: rows (a,b) -= q*(c,d)
        a, b = a - q * c, b - q * d
        if q:
            toks.append(('T', q))
        # now |a| < |c| (0<=a<|c| or so); left multiply by S^{-1} = [[0,1],[-1,0]]
        a, b, c, d = c, d, -a, -b
        toks.append(('S', 1))
    # now M' = [[a,b],[0,d]] with a*d = 1 -> +-T^{b*a}
    toks.append(('T', b * a))
    return toks


def check_st_word(M):
    toks = st_word(M)
    R = I2
    T = X
    for t, k in toks:
        if t == 'T':
            P = I2
            base = T if k >= 0 else inv(T)
            for _ in range(abs(k)):
                P = mul(P, base)
            R = mul(R, P)
        else:
            R = mul(R, S)
    return eq_pm(R, M)


# ---------- conjugacy-class canonical form of hyperbolic elements ----------
def _cf_period(a, b, c, d):
    """Attracting fixed point of z->(az+b)/(cz+d) (tr>2), regular CF tail period."""
    t = a + d
    D = t * t - 4
    sD = isqrt(D)
    assert sD * sD != D, 'square discriminant'
    # fixed points z = ((a-d) +- sqrt(D)) / (2c)
    import mpmath
    mpmath.mp.dps = 60
    best = None
    for s in (1, -1):
        z = (mpmath.mpf(a - d) + s * mpmath.sqrt(D)) / (2 * c)
        der = 1 / (c * z + d) ** 2
        if abs(der) < 1:
            best = s
    s = best
    # z = (P + sqrt(D'))/Q form with Q | D' - P^2 : use z = ((a-d)*s' + s sqrt(D))/(2c)
    P, Q, Dd = (a - d), 2 * c, D
    if s == -1:
        P, Q = -P, -Q  # z = (-(a-d) + sqrt D)/(-2c)
    # ensure Q | Dd - P^2
    if (Dd - P * P) % Q != 0:
        P, Q, Dd = P * abs(Q), Q * abs(Q), Dd * Q * Q
    seen = {}
    terms = []
    sq = mpmath.sqrt(Dd)
    while (P, Q) not in seen:
        seen[(P, Q)] = len(terms)
        x = (P + sq) / Q
        a0 = int(mpmath.floor(x))
        terms.append(a0)
        P = a0 * Q - P
        Q = (Dd - P * P) // Q
    start = seen[(P, Q)]
    per = terms[start:]
    if len(per) % 2:
        per = per + per
    return per, start % 2


def canon_class(A):
    """Canonical positive XY-word (min over rotations) of a hyperbolic class of PSL2(Z)."""
    if tr(A) < 0:
        A = neg(A)
    assert tr(A) > 2
    a, b, c, d = A
    if c == 0:
        raise ValueError
    per, par = _cf_period(a, b, c, d)
    w = ''.join(('X' if (i + par) % 2 == 0 else 'Y') * n for i, n in enumerate(per))
    w = primitive_root(w)
    n = len(w)
    return min((w + w)[i:i + n] for i in range(n))


def canon_word(w):
    n = len(w)
    return min((w + w)[i:i + n] for i in range(n))


def primitive_root(w):
    n = len(w)
    for p in range(1, n + 1):
        if n % p == 0 and w[:p] * (n // p) == w:
            return w[:p]
