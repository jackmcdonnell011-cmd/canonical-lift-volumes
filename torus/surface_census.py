"""Closed curves on S listed by word length in the free basis of pi_1(S) (surface_words.py).

usage: python -I surface_census.py <surface> <N> [out.csv]

For each unoriented closed curve (cyclic word u, |u| <= N): its modular code word w and m (the curve
is the lift of w with code word w^m), the coset cycle of that lift, and whether it fills S (exact test,
filling_exact.py).  Peripheral curves (around a puncture) are listed but cannot fill.
"""
import sys, os, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
surface, N = sys.argv[1], int(sys.argv[2])
out = sys.argv[3] if len(sys.argv) > 3 else None
sys.argv = [sys.argv[0], surface, os.devnull]
import filling_exact as fe
import psl2
from lifts import _cycles
from surface_words import Basis, enumerate_curves, cyclic_canon, inverse, letters

G = fe.G
B = Basis(G)


def mpow(M, k):
    P = psl2.I2
    for _ in range(k):
        P = psl2.mul(P, M)
    return P


def locate(u):
    """(w, m, c, orientation) with the curve u = lift of w through coset cycle containing c."""
    g = B.matrix(u)
    t = abs(psl2.tr(g))
    if t == 2:
        return None
    w = psl2.canon_class(g)
    Mw = psl2.xy_matrix(w)
    m = 1
    while abs(psl2.tr(mpow(Mw, m))) < t:
        m += 1
    assert abs(psl2.tr(mpow(Mw, m))) == t
    target, target_inv = cyclic_canon(u), cyclic_canon(inverse(u))
    for cy in _cycles(G.perm(Mw)):
        if len(cy) != m:
            continue
        c = cy[0]
        h = G.reps[c]
        dlt = psl2.mul(psl2.mul(h, mpow(Mw, m)), psl2.inv(h))
        k = cyclic_canon(B.word(dlt))
        if k == target:
            return (w, m, c, +1)
        if k == target_inv:
            return (w, m, c, -1)
    raise RuntimeError(f'no lift of {w}^{m} matches {letters(u)}')


if __name__ == '__main__':
    rows = []
    tot, fil, per = collections.Counter(), collections.Counter(), collections.Counter()
    for u in enumerate_curves(len(B.gens), N):
        L = len(u)
        tot[L] += 1
        loc = locate(u)
        if loc is None:
            per[L] += 1
            rows.append(dict(surface_word=letters(u), surface_word_length=L, peripheral=True))
            continue
        w, m, c, o = loc
        f = fe.fills(w, c, m)
        fil[L] += bool(f)
        rows.append(dict(surface_word=letters(u), surface_word_length=L, peripheral=False, modular_word=w,
                         m=m, code_length=m * len(w), coset=c, orientation=o, fills=f))
    print(f'{surface}: generators {len(B.gens)}; curves by surface word length (filling / peripheral / all)')
    for L in sorted(tot):
        fl = [r['code_length'] for r in rows if r['surface_word_length'] == L and r.get('fills')]
        rng = f'  code length of filling ones {min(fl)}..{max(fl)}' if fl else ''
        print(f'  {L}: {fil[L]} / {per[L]} / {tot[L]}{rng}')
    if out:
        cols = ['surface_word', 'surface_word_length', 'peripheral', 'modular_word', 'm', 'code_length',
                'coset', 'orientation', 'fills']
        with open(out, 'w', newline='') as f:
            wr = csv.DictWriter(f, fieldnames=cols)
            wr.writeheader()
            wr.writerows(rows)
