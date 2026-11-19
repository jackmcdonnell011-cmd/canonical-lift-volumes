"""Parent manifolds of the twist families X^kY, X^kY^2, X^kY^3 (Section 3, Table of parents).

For each family:
  1. Drill the shortest closed geodesic out of M_w for k = 10.  This gives a three-cusped manifold N
     (cusp 0 = Lorenz knot K, cusp 1 = trefoil T, cusp 2 = the drilled curve c).  Drilling the members
     with k = 7 and k = 14 is checked to give an isometric N.
  2. Filling K and T by their meridians gives S^3 minus c with fundamental group Z, so c is an unknot.
     Its Seifert longitude lambda_c is the slope whose filling also has fundamental group Z, and H_1 of
     the 0-surgery on c shows that c links K and T once each.
  3. The fillings of c along mu_c + n lambda_c (the 1/n surgeries on the unknot c) are compared with the
     members k = 1..20 of the family, by volume and SnapPy's is_isometric_to.

usage:  python twist_parents.py [out_dir]     writes twist_parents.csv (one row per filling that is a family
        member) and twist_parents_summary.txt
"""
import csv
import os
import sys

import snappy
import spherogram

from dt_codes import canonical, get_dt

FAMILIES = {
    'X^kY': lambda k: 'X' * k + 'Y',
    'X^kY^2': lambda k: 'X' * k + 'YY',
    'X^kY^3': lambda k: 'X' * k + 'YYY',
}
K_DRILL = (10, 7, 14)
K_MAX = 20


def positive(M, tries=20):
    """re-triangulate until SnapPy finds a solution with all tetrahedra positively oriented"""
    for _ in range(tries):
        if M.solution_type() == 'all tetrahedra positively oriented':
            return M
        M.randomize()
    return M


def self_crossings(L):
    comp = {}
    for i, c in enumerate(L.link_components):
        for cs in c:
            comp.setdefault(cs.crossing, set()).add(i)
    cnt = [0] * len(L.link_components)
    for C in L.crossings:
        s = comp[C]
        if len(s) == 1:
            cnt[next(iter(s))] += 1
    return cnt


def manifold_of(w):
    """M_w with cusp 0 = Lorenz knot of w and cusp 1 = trefoil"""
    L = spherogram.Link('DT: ' + get_dt(canonical(w)))
    assert self_crossings(L)[1] == 3, 'component 1 should be the trefoil'
    return positive(L.exterior())


def pi1_is_Z(M):
    G = M.fundamental_group()
    return G.num_generators() == 1 and len(G.relators()) == 0


def main(out_dir):
    out, rows = [], []
    for name, f in FAMILIES.items():
        out.append(f'\n===== {name} =====')
        fam = {k: manifold_of(f(k)) for k in range(1, K_MAX + 1)}
        vol = {k: float(fam[k].volume()) for k in fam}

        drilled = {k0: positive(manifold_of(f(k0)).drill(0)) for k0 in K_DRILL}
        N = drilled[K_DRILL[0]]
        out.append(f'N (drilled from k={K_DRILL[0]}): {N.num_cusps()} cusps, volume {float(N.volume()):.10f}, '
                   f'{N.solution_type()}, census identification {N.identify()}')
        for k0 in K_DRILL[1:]:
            out.append(f'  drilling k={k0} gives an isometric manifold: {N.is_isometric_to(drilled[k0])}')

        S = N.copy()
        S.dehn_fill([(1, 0), (1, 0), (0, 0)])
        out.append(f'  N(mu_K, mu_T) = S^3 minus c has fundamental group Z: {pi1_is_Z(S)}')
        lam = None
        for p in range(-40, 41):
            T = N.copy()
            T.dehn_fill([(1, 0), (1, 0), (p, 1)])
            if pi1_is_Z(T) and str(T.homology()) == 'Z':
                lam = (p, 1)
                break
        out.append(f'  Seifert longitude of c: {lam}')
        if lam is None:
            continue
        for nm, keep in (('K', 0), ('T', 1)):
            T = N.copy()
            fills = [(1, 0), (1, 0), lam]
            fills[keep] = (0, 0)
            T.dehn_fill(fills)
            out.append(f'  H_1(S^3_0(c) minus {nm}) = {T.homology()}   (Z means c links {nm} once)')

        hits = {}
        for n in range(-25, 26):
            F = N.copy()
            F.dehn_fill((1 + n * lam[0], n * lam[1]), 2)
            F = positive(F)
            v = float(F.volume())
            iso = [k for k in fam if abs(vol[k] - v) < 1e-7 and F.is_isometric_to(fam[k])]
            if iso:
                hits[n] = iso[0]
                rows.append(dict(family=name, n=n, slope=f'{(1 + n * lam[0], n * lam[1])}', k=iso[0],
                                 volume=f'{v:.10f}', parent_volume=f'{float(N.volume()):.10f}'))
        missing = [k for k in range(1, K_MAX + 1) if k not in hits.values()]
        out.append('  fillings that are family members: ' + ', '.join(f'n={n} -> k={k}' for n, k in sorted(hits.items())))
        out.append(f'  members k=1..{K_MAX} reproduced: {K_MAX - len(missing)}; missing: {missing}')
        out.append(f'  volumes: k=1 {vol[1]:.4f}, k=10 {vol[10]:.4f}, k=20 {vol[20]:.4f}, vol(N) {float(N.volume()):.4f}; '
                   f'increasing in k: {all(vol[k] < vol[k + 1] for k in range(1, K_MAX))}')

    txt = '\n'.join(out)
    print(txt)
    os.makedirs(out_dir, exist_ok=True)
    open(os.path.join(out_dir, 'twist_parents_summary.txt'), 'w').write(txt.lstrip() + '\n')
    with open(os.path.join(out_dir, 'twist_parents.csv'), 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
