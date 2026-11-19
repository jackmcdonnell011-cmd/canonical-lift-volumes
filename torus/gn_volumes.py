"""volumes(Gamma, L): volumes of UT(S) - (canonical lift of kappa) for every oriented
closed geodesic kappa of length <= L on S = H^2/Gamma, Gamma <= PSL2(Z) torsion-free.

Usage:  python -I gn_volumes.py torus11 <L> [out.csv]
  torus11: the once-punctured torus, S = H^2/Gamma' with Gamma' the commutator subgroup (index 6)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import csv
import math
import time
import psl2
from lifts import ModularLink, word_length
import subgroups as sg

SURFACES = {
    'torus11': lambda: sg.Commutator(),
}


def modular_words(L):
    """All primitive cyclic words in X,Y (both letters, canonical rotation) with
    hyperbolic length <= L on the modular surface (one per oriented closed geodesic)."""
    tmax = 2 * math.cosh(L / 2)
    out = []

    def dfs(w, A):
        if 'X' in w and 'Y' in w and w == psl2.canon_word(w) and psl2.primitive_root(w) == w:
            out.append(w)
        for ch, B in (('X', psl2.X), ('Y', psl2.Y)):
            C = psl2.mul(A, B)
            # trace only grows when letters are added (all matrices >= I entrywise);
            # a word must eventually contain Y, so bound by tr(C*Y) while Y is absent
            Cb = C if 'Y' in w + ch else psl2.mul(C, psl2.Y)
            if psl2.tr(Cb) <= tmax + 1e-9:
                dfs(w + ch, C)
    dfs('X', psl2.X)  # canonical rotations of words containing X start with X
    return sorted(set(out), key=lambda w: (word_length(w), w))


def hyp_volume(N, tries=20):
    """Volume of a (possibly Dehn-filled) manifold, retriangulating until the
    solution has all tetrahedra positively oriented. Returns (vol, solution_type)."""
    F = N.filled_triangulation() if any(not c['complete?'] for c in N.cusp_info()) else N.copy()
    st = F.solution_type()
    for _ in range(tries):
        if st == 'all tetrahedra positively oriented':
            break
        F.randomize()
        st = F.solution_type()
    return float(F.volume()), st


def volumes(Gamma, L, verbose=True, log=None):
    sig = Gamma.signature()
    assert sig['e2'] == 0 and sig['e3'] == 0, 'Gamma has torsion'
    rows = []
    words = modular_words(L)
    if verbose:
        print(f'# {Gamma.name}: genus {sig["genus"]}, cusps {sig["cusps"]}, index {Gamma.d}; '
              f'{len(words)} modular words with length <= {L}', flush=True)
    for w in words:
        lw = word_length(w)
        # cheap prefilter: lifts of w are the cycles of w acting on the cosets
        p = Gamma.perm(psl2.xy_matrix(w))
        seen, mins = set(), []
        for x in range(Gamma.d):
            if x in seen:
                continue
            j, m = x, 0
            while j not in seen:
                seen.add(j)
                j = p[j]
                m += 1
            mins.append(m)
        if min(mins) * lw > L + 1e-9:
            continue
        t0 = time.time()
        ml = ModularLink(w)
        cd = ml.cover(Gamma)
        fullvol = float(cd.C.volume())
        for (i, m, N) in cd.single_lifts():
            if m * lw > L + 1e-9:
                continue
            v, st = hyp_volume(N)
            row = dict(surface=Gamma.name, g=sig['genus'], n=sig['cusps'], word=w, m=m, lift_cusp=i,
                       geodesic=f'({w})^{m}' if m > 1 else w, length=m * lw, volume=v, solution=st,
                       modular_volume=ml.vol, full_preimage_ratio=fullvol / ml.vol)
            rows.append(row)
            if verbose:
                print(f'{w:>14} m={m} len={m*lw:8.4f} vol={v:10.6f} {st}', flush=True)
        if verbose:
            print(f'   [{w}: {len(cd.Kcusps)} lifts, full preimage vol/modular vol = {fullvol/ml.vol:.9f}, '
                  f'{time.time()-t0:.1f}s]', flush=True)
    return rows


if __name__ == '__main__':
    name, L = sys.argv[1], float(sys.argv[2])
    rows = volumes(SURFACES[name](), L)
    if len(sys.argv) > 3:
        with open(sys.argv[3], 'w', newline='') as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            wr.writeheader()
            wr.writerows(rows)
