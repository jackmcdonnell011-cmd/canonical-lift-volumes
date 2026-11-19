"""Checked pipeline: single canonical-lift volumes on S = H^2/Gamma, with checks 1-5,
exact lift <-> Gamma-element bookkeeping, symmetry/isometry de-duplication.

Usage: python -I pipeline.py <surface> <L>
Writes  results_<surface>_L<L>.csv  (one row per geodesic of S up to symmetry)
        lifts_<surface>_L<L>.csv    (one row per lifted component, before de-duplication)
        summary_<surface>_L<L>.txt
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import csv
import math
import time
import psl2
import subgroups as sg
from lifts import ModularLink, word_length, _act_word, _cycles, _inv
from gn_volumes import SURFACES, modular_words
import snappy

POS = 'all tetrahedra positively oriented'


# ---------------------------------------------------------------- PSL2(Z) helpers
def mpow(A, k):
    R = psl2.I2
    B = A if k >= 0 else psl2.inv(A)
    for _ in range(abs(k)):
        R = psl2.mul(R, B)
    return R


def conjugator(A, w):
    """g with g^{-1} A g = +-M_w (w a canonical positive word), else None."""
    Mw = psl2.xy_matrix(w)
    if psl2.tr(A) < 0:
        A = psl2.neg(A)
    a, b, c, d = A
    per, par = psl2._cf_period(a, b, c, d)
    # regenerate the full CF expansion (pre-period + periods) by recomputing terms
    import mpmath
    mpmath.mp.dps = 80
    D = (a + d) ** 2 - 4
    best = None
    for s in (1, -1):
        z = (mpmath.mpf(a - d) + s * mpmath.sqrt(D)) / (2 * c)
        if abs(1 / (c * z + d) ** 2) < 1:
            best = z
    x = best
    terms = []
    for _ in range(4 * len(w) + 60):
        a0 = int(mpmath.floor(x))
        terms.append(a0)
        fr = x - a0
        if fr == 0:
            break
        x = 1 / fr
    # candidate prefixes: block X^{a0}, then letters one at a time
    P = mpow(psl2.X, terms[0])
    letters = []
    for i, t in enumerate(terms[1:], start=1):
        letters += [psl2.Y if i % 2 else psl2.X] * t
    cands = [P]
    for L in letters[:6 * len(w) + 40]:
        P = psl2.mul(P, L)
        cands.append(P)
    for P in cands:
        for Q in (P, psl2.mul(P, psl2.S)):
            B = psl2.mul(psl2.mul(psl2.inv(Q), A), Q)
            if psl2.eq_pm(B, Mw):
                return Q
    # try rotations of w explicitly
    n = len(w)
    for P in cands:
        B = psl2.mul(psl2.mul(psl2.inv(P), A), P)
        for r in range(n):
            pre = psl2.xy_matrix(w[:r])
            Q = psl2.mul(P, psl2.inv(pre)) if False else None
        # rotation: M_{w[r:]+w[:r]} = pre^{-1} M_w pre
        for r in range(n):
            rot = w[r:] + w[:r]
            if psl2.eq_pm(B, psl2.xy_matrix(rot)):
                pre = psl2.xy_matrix(w[:r])
                Q = psl2.mul(P, psl2.inv(pre))
                assert psl2.eq_pm(psl2.mul(psl2.mul(psl2.inv(Q), A), Q), Mw)
                return Q
    return None


def inverse_word(w):
    """canonical word of the reversed geodesic: M_w^{-1} ~ M_{rev(swap(w))}"""
    sw = w.translate(str.maketrans('XY', 'YX'))
    return psl2.canon_word(sw[::-1])


def swap_word(w):
    return psl2.canon_word(w.translate(str.maketrans('XY', 'YX')))


# ---------------------------------------------------------------- Gamma structure
def normaliser_autos(G):
    """Automorphisms of the transitive PSL2(Z)-set (= N(Gamma)/Gamma acting on cosets)."""
    d = G.d
    out = []
    for y in range(d):
        phi = {0: y}
        stack = [0]
        ok = True
        while stack and ok:
            x = stack.pop()
            for p in (G.pS, G.pT):
                a, b = p[x], p[phi[x]]
                if a in phi:
                    if phi[a] != b:
                        ok = False
                        break
                else:
                    phi[a] = b
                    stack.append(a)
        if ok and len(set(phi.values())) == d:
            out.append([phi[i] for i in range(d)])
    return out


def J_map(G):
    """If J Gamma J = Gamma exactly, return the induced map on cosets x -> coset of J h_x J."""
    for x in range(G.d):
        for g in (psl2.S, psl2.X):
            y = G.act(x, g)
            gam = psl2.mul(psl2.mul(G.reps[x], g), psl2.inv(G.reps[y]))  # Schreier generator
            if G.act(0, psl2.conjJ(gam)) != 0:
                return None
    return [G.act(0, psl2.conjJ(G.reps[x])) for x in range(G.d)]


def trefoil_cover(G):
    """Cover of S^3 - trefoil for Gamma, built from SnapPy's 3_1 with a -> S, b -> ST (no K)."""
    T = snappy.Manifold('3_1')
    FG = T.fundamental_group()
    assert FG.relators() == ['aabbb'], FG.relators()
    U = psl2.mul(psl2.S, psl2.X)
    gm = {'a': psl2.S, 'b': U}
    assert psl2.eq_pm(psl2.word_matrix('aabbb', gm), psl2.I2)
    perms = [G.perm(gm['a']), G.perm(gm['b'])]
    try:
        return T.cover(perms)
    except Exception:
        return T.cover([_inv(p) for p in perms])


def _abelian(rows, n):
    """(betti number, torsion) of Z^n modulo the integer vectors in rows.  Generators with a coefficient +-1 in
    some row are eliminated first (sparse, exact); the small remainder goes to PARI's Smith normal form."""
    from cypari import pari
    R = [{j: v for j, v in enumerate(r) if v} for r in rows]
    R = [r for r in R if r]
    alive = set(range(n))
    while True:
        best = None                               # unit pivot in the shortest row, to limit fill-in
        for i, r in enumerate(R):
            if best is None or len(r) < len(R[best[0]]):
                j = next((j for j, v in r.items() if v in (1, -1)), None)
                if j is not None:
                    best = (i, j)
        if best is None:
            break
        i, j = best
        r = R.pop(i)
        v = r[j]                                  # v = +-1, so r says x_j = -v * (rest of r)
        for s in R:
            c = s.get(j)
            if c:
                f = c * v
                for k, a in r.items():
                    x = s.get(k, 0) - f * a
                    if x:
                        s[k] = x
                    else:
                        s.pop(k, None)
        alive.discard(j)
        R = [s for s in R if s]
    keep = sorted(alive)
    m = len(keep)
    if not R:
        return (m, ())
    N = max(len(R), m)
    flat = []
    for i in range(N):
        flat.extend([R[i].get(j, 0) for j in keep] + [0] * (N - m) if i < len(R) else [0] * N)
    d = [abs(int(x)) for x in pari.matsnf(pari.matrix(N, N, flat))]
    return (sum(1 for x in d if x == 0) - (N - m), tuple(sorted(x for x in d if x > 1)))


def _homs_mod2(rels, n):
    """basis (bitmasks over the generators) of the homomorphisms from <gens | rels> onto Z/2"""
    pivots = {}                                   # pivot bit -> reduced relator parity row
    for r in rels:
        v = 0
        for x in r:
            v ^= 1 << (abs(x) - 1)
        for p, row in pivots.items():
            if v >> p & 1:
                v ^= row
        if v:
            p = v.bit_length() - 1
            for q in list(pivots):
                if pivots[q] >> p & 1:
                    pivots[q] ^= v
            pivots[p] = v
    basis = []
    for f in range(n):                            # one solution per free generator
        if f in pivots:
            continue
        phi = 1 << f
        for p, row in pivots.items():
            if row >> f & 1:
                phi |= 1 << p
        basis.append(phi)
    return basis


def _double_cover_h1(rels, n, phi):
    """H1 of the index-2 subgroup ker(phi), by Reidemeister-Schreier with transversal {1, t}"""
    t = next(i for i in range(n) if phi >> i & 1)
    col = {}
    for g in range(n):
        for c in (0, 1):
            if (g, c) != (t, 0):                  # y(t, 0) = t t^-1 is trivial
                col[(g, c)] = len(col)
    rows = []
    for r in rels:
        for c in (0, 1):
            v = [0] * len(col)
            for x in r:
                g = abs(x) - 1
                if x > 0:
                    if (g, c) in col:
                        v[col[(g, c)]] += 1
                    c ^= phi >> g & 1
                else:
                    c ^= phi >> g & 1
                    if (g, c) in col:
                        v[col[(g, c)]] -= 1
            rows.append(v)
    return _abelian(rows, len(col))


def cover_fingerprint(N):
    """(unfilled cusps, H1, sorted H1 of every double cover); N may have Dehn-filled cusps.
    Computed from the unsimplified presentation of pi_1 (which includes the filling relations),
    by Reidemeister-Schreier and Smith normal forms.  SnapPy's own homology() / covers(2) simplify pi_1 first,
    which can fail to finish on large filled covers."""
    G = N.fundamental_group(simplify_presentation=False)
    n = len(G.generators())
    rels = G.relators(as_int_list=True)
    H = _abelian([[sum((x > 0) - (x < 0) for x in r if abs(x) == g + 1) for g in range(n)] for r in rels], n)
    basis = _homs_mod2(rels, n)
    phis = set()
    for k in range(1, 2 ** len(basis)):
        phi = 0
        for i, b in enumerate(basis):
            if k >> i & 1:
                phi ^= b
        phis.add(phi)
    c2 = sorted(_double_cover_h1(rels, n, phi) for phi in phis)
    unfilled = sum(1 for f in N.cusp_info('filling') if tuple(f) == (0, 0))
    return (unfilled, H, tuple(c2))


# ---------------------------------------------------------------- per word
def hyp(N, tries=20):
    F = N.filled_triangulation()
    st = F.solution_type()
    for _ in range(tries):
        if st == POS:
            break
        F.randomize()
        st = F.solution_type()
    return F, float(F.volume()), st


_BASES = {}
def _basis(G):
    import surface_words as sw
    if id(G) not in _BASES:
        _BASES[id(G)] = sw.Basis(G)
    return _BASES[id(G)]


def process_word(w, G, L, Dfp, sig, cache_dir=None):
    res = dict(word=w, checks={}, rows=[], notes=[])
    ml = ModularLink(w, cache_dir=cache_dir)
    gm, Gp = ml.gm, ml.G
    per = Gp.peripheral_curves()
    muK, lamK = per[ml.cuspK]
    muT, lamT = per[ml.cuspT]
    # ---- check 1: rho kills mu_K and the fibre; rho irreducible; mu_T -> X^{+-1}
    c1 = psl2.eq_pm(psl2.word_matrix(muK, gm), psl2.I2)
    fib = None
    for sgn in (6, -6):
        f = (muT * sgn if sgn > 0 else muT.swapcase()[::-1] * 6) + lamT
        if psl2.eq_pm(psl2.word_matrix(f, gm), psl2.I2):
            fib = f'mu_T^{sgn} lambda_T'
    c1 = c1 and fib is not None
    mT = psl2.word_matrix(muT, gm)
    c1 = c1 and (psl2.eq_pm(mT, psl2.X) or psl2.eq_pm(mT, psl2.inv(psl2.X)))
    irr = any(psl2.tr(psl2.mul(psl2.mul(A, B), psl2.mul(psl2.inv(A), psl2.inv(B)))) not in (2, -2)
              for A in gm.values() for B in gm.values())
    c1 = c1 and irr
    for r in Gp.relators():
        c1 = c1 and psl2.eq_pm(psl2.word_matrix(r, gm), psl2.I2)
    res['checks']['1_rho'] = c1
    res['fibre'] = fib
    # ---- check 3: core of K -> conjugate of +M_w, exact
    lam = psl2.word_matrix(lamK, gm)
    lamo = lam if ml.Korient == 1 else psl2.inv(lam)
    g0 = conjugator(lamo, w)
    Mw = psl2.xy_matrix(w)
    c3 = g0 is not None and abs(psl2.tr(lamo)) == psl2.tr(Mw)
    res['checks']['3_core'] = c3
    res['Korient'] = ml.Korient
    res['twisted'] = ml.twisted
    res['solver'] = ml.solver
    # ---- cover
    cd = ml.cover(G)
    C = cd.C
    d = G.d
    # tetrahedron correspondence -> cusps over K / T and sheets
    bidx = ml.M._get_cusp_indices_and_peripheral_curve_data()[0]
    cidx = C._get_cusp_indices_and_peripheral_curve_data()[0]
    nb = len(bidx)
    over, sheets = {}, {}
    for i in range(len(cidx)):
        for v in range(4):
            c = cidx[i][v]
            over.setdefault(c, bidx[i // d][v])
            sheets.setdefault(c, set()).add(i % d)
    Tcusps = [c for c in over if over[c] == ml.cuspT]
    Kcusps = [c for c in over if over[c] == ml.cuspK]
    # ---- check 2: topology of S from the cover
    n_cov = len(Tcusps)
    chi = -d // 6
    g2 = 2 - n_cov - chi
    c2 = (len(cidx) == nb * d and abs(float(C.volume()) / ml.vol - d) < 1e-8 and g2 % 2 == 0
          and n_cov == sig['cusps'] and g2 // 2 == sig['genus'])
    res['checks']['2_topology'] = c2
    res['g_cover'], res['n_cover'] = g2 // 2, n_cov
    # ---- check 4: degrees of K-lift cusps (from the cover's tetrahedra) must equal the
    # cycle lengths of M_w on Gamma\\PSL2(Z) (double cosets, computed independently);
    # each cycle gives a Gamma element h M_w^m h^-1 whose trace must be tr(M_w^m).
    base_cnt = sum(1 for t in range(nb) for v in range(4) if bidx[t][v] == ml.cuspK)
    cnt = {}
    for i in range(len(cidx)):
        for v in range(4):
            if bidx[i // d][v] == ml.cuspK:
                cnt[cidx[i][v]] = cnt.get(cidx[i][v], 0) + 1
    pw = G.perm(Mw)
    cyc_w = _cycles(pw)
    degs_indep = sorted(len(c) for c in cyc_w)
    kc = {i: (m, sl) for (i, m, sl) in cd.Kcusps}
    degs_cover = sorted(cnt[c] // base_cnt for c in Kcusps)
    ok4 = degs_cover == degs_indep and all(cnt[c] // base_cnt == kc[c][0] for c in Kcusps)
    gam_by_m = {}
    for cy in cyc_w:
        m = len(cy)
        h = G.reps[min(cy)]
        gam = psl2.mul(psl2.mul(h, mpow(Mw, m)), psl2.inv(h))
        ok4 = ok4 and G.act(0, gam) == 0 and abs(psl2.tr(gam)) == abs(psl2.tr(mpow(Mw, m)))
        gam_by_m.setdefault(m, []).append(gam)
    lifts = []
    for c in sorted(Kcusps):
        m = cnt[c] // base_cnt
        gam = gam_by_m[m][0]
        lifts.append(dict(cusp=c, m=m, gamma=gam, trace=abs(psl2.tr(gam)),
                          gamma_pinned=len(gam_by_m[m]) == 1, length=m * ml.length, slope=kc[c][1]))
    res['checks']['4_lifts'] = ok4
    res['degrees'] = degs_indep
    # ---- check 5: fill every K-lift with its meridian -> must be the trefoil cover for Gamma
    A = cd.all_filled()
    # both tests use the Dehn-filled cover A itself rather than a filled triangulation of it (the same manifold),
    # which avoids searching for a hyperbolic structure on UT(S), a Seifert fibred space.
    fp = cover_fingerprint(A)
    nonhyp = A.solution_type() != POS or float(A.volume()) < 1e-6
    # meridians are lifts of mu_K (exact integrality check) and agree with shape matching
    c5 = (fp == Dfp) and nonhyp and cd.meridians_ok
    res['checks']['5_meridians'] = c5
    res['meridian_method'] = 'lifted mu_K (checked against cusp shapes)'
    # ---- which curve on S each K-lift is: the Gamma element of its core, as an unoriented cyclic word
    # in the free basis of pi_1(S) (surface_words.py).  Sheet x of base tetrahedron t lies in cover
    # tetrahedron t*d + x, but each base tetrahedron is reached from the basepoint by its own path, so
    # only at some (t, vertex) do the sheets split into K-cusps exactly like the cycles of lambda_K.
    # Read the curves off every such (t, vertex), for lambda_K from both SnapPy's simplified and
    # unsimplified presentations (the unsimplified one is tied to the tetrahedra); keep the answer only
    # if they all agree, else leave it unknown (None).
    import surface_words as sw
    B = _basis(G)
    G0 = ml.M.fundamental_group(simplify_presentation=False)
    om = dict(zip(G0.generators(), [psl2.word_matrix(x, gm) for x in Gp.original_generators()]))
    lams = [lam, psl2.word_matrix(G0.peripheral_curves()[ml.cuspK][1], om)]
    Kset = set(Kcusps)
    readings = set()
    for lam_x in lams:
        lam_cyc = {frozenset(c) for c in _cycles(G.perm(lam_x))}
        for t in range(nb):
            for v in range(4):
                if bidx[t][v] != ml.cuspK:
                    continue
                part = {}
                for x in range(d):
                    part.setdefault(cidx[t * d + x][v], set()).add(x)
                if {frozenset(s) for s in part.values()} != lam_cyc or not set(part) <= Kset:
                    continue
                reading = []
                for c, s in part.items():
                    x = min(s)
                    core = psl2.mul(psl2.mul(G.reps[x], mpow(lam_x, len(s))), psl2.inv(G.reps[x]))
                    if G.act(0, core) != 0:
                        reading = None
                        break
                    reading.append((c, sw.letters(sw.unoriented_canon(B.word(core)))))
                if reading:
                    readings.add(tuple(sorted(reading)))
    res['surface_word_readings'] = len(readings)
    # ---- single lifts
    for lf in lifts:
        if lf['m'] * len(w) > L:
            continue
        N = C.copy()
        # all the other K-lifts in one dehn_fill call: SnapPea re-solves the whole cover after every call
        fill = [(0, 0)] * C.num_cusps()
        for (j, mj, sl) in cd.Kcusps:
            if j != lf['cusp']:
                fill[j] = sl
        N.dehn_fill(fill)
        F, v, st = hyp(N)
        remaining = sorted(j for j in range(C.num_cusps()) if j == lf['cusp'] or j in Tcusps)
        lf.update(volume=v, solution=st, F=F, Kpos=remaining.index(lf['cusp']))
        res['rows'].append(lf)
    # If the readings disagree (a symmetry of the cover can make a wrong (t, vertex) split the sheets the
    # same way), keep only readings where every hyperbolic lift is a curve that fills S (exact test,
    # filling_exact.py): a curve that does not fill never has a hyperbolic complement.
    if len(readings) > 1:
        import filling_exact as fe
        if fe.G is not G:
            fe.setup(G)
        fills = {}
        for cy in cyc_w:
            m = len(cy)
            core = psl2.mul(psl2.mul(G.reps[cy[0]], mpow(Mw, m)), psl2.inv(G.reps[cy[0]]))
            fills[sw.letters(sw.unoriented_canon(B.word(core)))] = fe.fills(w, cy[0], m)
        hypc = {lf['cusp'] for lf in res['rows'] if lf['solution'] == POS}
        readings = {rd for rd in readings if all(fills.get(wd) for c, wd in rd if c in hypc)}
        res['surface_word_readings_after_filling'] = len(readings)
        # still several, but all give every curve the same volume: the curve -> volume data agree
        vol = {lf['cusp']: (round(lf['volume'], 6) if lf['solution'] == POS else None) for lf in res['rows']}
        # (non-hyperbolic cusps have vol None, so sort on (word, None-flag, volume))
        key = lambda x: (x[0], x[1] is None, x[1] or 0.0)
        if len(readings) > 1 and len({tuple(sorted(((wd, vol.get(c)) for c, wd in rd), key=key))
                                      for rd in readings}) == 1:
            readings = {min(readings)}
            res['surface_word_note'] = 'readings differ only between curves of equal volume'
    # exact labels (lift_labels.py) when its corner-by-corner self-check passes; the readings above are
    # kept only as a fallback.
    from lift_labels import exact_labels
    ex = exact_labels(ml, G, B, C)
    if ex is not None:
        cusp_word = ex
        res['label_method'] = 'exact'
    else:
        cusp_word = dict(next(iter(readings))) if len(readings) == 1 else {}
        res['label_method'] = 'corner readings'
    # Gamma element of each labelled lift from its own coset cycle (gam_by_m[m][0] above is only one of the
    # equal-m components); a reciprocal w has a curve and its reverse with the same unoriented label
    gam_by_label = {}
    for cy in cyc_w:
        m = len(cy)
        gam = psl2.mul(psl2.mul(G.reps[cy[0]], mpow(Mw, m)), psl2.inv(G.reps[cy[0]]))
        gam_by_label.setdefault((sw.letters(sw.unoriented_canon(B.word(gam))), m), gam)
    for lf in lifts:
        lf['surface_word'] = cusp_word.get(lf['cusp'])
        gam = gam_by_label.get((lf['surface_word'], lf['m']))
        if gam is not None:
            lf.update(gamma=gam, trace=abs(psl2.tr(gam)), gamma_pinned=True)
    res['ml_vol'] = ml.vol
    res['ml_len'] = ml.length
    return res


def cusp_iso(F1, k1, F2, k2):
    try:
        isos = F1.is_isometric_to(F2, return_isometries=True)
    except Exception as e:
        return None
    for iso in isos:
        if iso.cusp_images()[k1] == k2:
            return True
    return False


def run(surface, L):
    G = SURFACES[surface]()
    sig = G.signature()
    autos = normaliser_autos(G)
    normal = len(autos) == G.d
    Jm = J_map(G)
    D = trefoil_cover(G)
    Dfp = cover_fingerprint(D)
    log = []
    P = lambda *a: (print(*a, flush=True), log.append(' '.join(str(x) for x in a)))
    P(f'# {surface}: {G.name}  index {G.d}  genus {sig["genus"]} cusps {sig["cusps"]} widths {sig["widths"]}; '
      f'normal in PSL2(Z): {normal}; J-invariant (exactly): {Jm is not None}; |N(Gamma)/Gamma| = {len(autos)}')
    P(f'# direct trefoil cover for Gamma: cusps {Dfp[0]}, H1 {Dfp[1]}')
    words = []
    for n in range(2, int(L) + 1):
        for k in range(2 ** n):
            w = ''.join('XY'[(k >> i) & 1] for i in range(n))
            if 'X' in w and 'Y' in w and w == psl2.canon_word(w) and psl2.primitive_root(w) == w:
                cyc = _cycles(G.perm(psl2.xy_matrix(w)))
                if min(len(c) for c in cyc) * len(w) <= L:
                    words.append(w)
    P(f'# {len(words)} modular words (|w| <= {L:g}) carry components with m*|w| <= {L:g}')
    results = {}
    lift_rows = []
    for w in words:
        t0 = time.time()
        try:
            r = process_word(w, G, L, Dfp, sig)
        except Exception as e:
            P(f'{w}: FAILED ({type(e).__name__}: {e})')
            results[w] = None
            continue
        results[w] = r
        ck = ''.join('P' if v else 'F' for v in r['checks'].values())
        P(f'{w}: checks(1-rho,3-core,2-topol,4-lifts,5-merid)={ck} degrees={r["degrees"]} '
          f'g,n from cover=({r["g_cover"]},{r["n_cover"]}) solver={r["solver"]} [{time.time()-t0:.1f}s]')
        for lf in r['rows']:
            P(f'    lift cusp {lf["cusp"]} m={lf["m"]} len={lf["length"]:.6f} trace={lf["trace"]} '
              f'vol={lf["volume"]:.10f} {lf["solution"][:22]}')
            lift_rows.append((w, lf))
    # ------------------------------------------------ de-duplication & isometry checks
    out_rows, suspects = [], []
    for w in words:
        r = results.get(w)
        if r is None:
            continue
        ck = ' '.join(f'{k}={"pass" if v else "FAIL"}' for k, v in r['checks'].items())
        allok = all(r['checks'].values())
        rows = r['rows']
        if not rows:
            continue
        if normal:
            ms = set(r['degrees'])
            status = []
            m = rows[0]['m']
            if len(ms) != 1 or len(r['degrees']) != G.d // m or len(rows) != G.d // m:
                status.append(f'SUSPECT: degrees {r["degrees"]} / {len(rows)} lifts, expected index/m equal lifts')
            hypw = all(lf['solution'] == POS for lf in rows)
            vols = [lf['volume'] for lf in rows]
            if hypw:
                if max(vols) - min(vols) > 1e-8:
                    status.append('SUSPECT: lifts of one word have different volumes')
                for lf in rows[1:]:
                    if cusp_iso(rows[0]['F'], rows[0]['Kpos'], lf['F'], lf['Kpos']) is not True:
                        status.append(f'SUSPECT: lift cusp {lf["cusp"]} not cusp-isometric to cusp {rows[0]["cusp"]}')
            groups = [(rows[0], len(rows), status, hypw)]
        else:
            groups = [(lf, 1, [], lf['solution'] == POS) for lf in rows]
        for lf, mult, status, hypc in groups:
            if not allok:
                status = status + ['EXCLUDED: a check failed']
            row = dict(word=w, modular_word_length=len(w), m=lf['m'], surface_word_length=lf['m'] * len(w),
                       hyperbolic_length_on_S=round(lf['length'], 10),
                       volume=round(lf['volume'], 10) if hypc else '',
                       multiplicity=mult, gamma_trace=lf['trace'], checks=ck,
                       gamma_element=str(lf['gamma']) + ('' if (lf['gamma_pinned'] or normal) else ' (one of the equal-m components)'),
                       hyperbolic=hypc, solution_type=lf['solution'],
                       status='; '.join(status) if status else 'ok', _lf=lf)
            if status:
                suspects.append((w, lf['cusp'], status))
            out_rows.append(row)
    # equal volumes between different rows -> cusp-respecting isometry test + relation label
    cross = []
    hr = sorted([r for r in out_rows if r['hyperbolic']], key=lambda r: r['_lf']['volume'])
    for i in range(len(hr)):
        for j in range(i + 1, len(hr)):
            A, B = hr[i], hr[j]
            if abs(A['_lf']['volume'] - B['_lf']['volume']) > 1e-8:
                break
            iso = cusp_iso(A['_lf']['F'], A['_lf']['Kpos'], B['_lf']['F'], B['_lf']['Kpos'])
            wa, wb = A['word'], B['word']
            if wa == wb:
                rel = 'same word, different component'
            elif inverse_word(wa) == wb:
                rel = 'inversion (kappa vs kappa^-1)'
            elif swap_word(wa) == wb:
                rel = 'J (mirror)' + ('' if Jm is not None else ', Gamma NOT J-invariant')
            elif psl2.canon_word(wa[::-1]) == wb:
                rel = 'J o inversion' + ('' if Jm is not None else ', Gamma NOT J-invariant')
            else:
                rel = 'no known symmetry'
            cross.append((wa, A['_lf']['cusp'], wb, B['_lf']['cusp'], round(A['_lf']['volume'], 9), iso, rel))
            for R, other in ((A, wb), (B, wa)):
                t = f'equal vol to {other}: cusp-isometry={iso}, relation={rel}'
                if iso is not True:
                    t = 'SUSPECT ' + t
                    suspects.append((R['word'], R['_lf']['cusp'], [t]))
                R['status'] = t if R['status'] == 'ok' else R['status'] + '; ' + t
    for r in out_rows:
        del r['_lf']
    out_rows.sort(key=lambda r: (r['surface_word_length'], r['hyperbolic_length_on_S'], r['word']))
    fn = f'results_{surface}_L{L:g}.csv'
    if out_rows:
        with open(fn, 'w', newline='') as f:
            wr = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
            wr.writeheader()
            wr.writerows(out_rows)
    clean = [r for r in out_rows if r['hyperbolic'] and 'SUSPECT' not in r['status'] and 'EXCLUDED' not in r['status']]
    distinct = len({round(r['volume'], 7) for r in clean})
    okw = [w for w in results if results[w]]
    checksum = {k: all(results[w]['checks'][k] for w in okw) for k in
                ['1_rho', '3_core', '2_topology', '4_lifts', '5_meridians']}
    P(f'# CHECKS over all processed words: {checksum}')
    P(f'# words processed {len(okw)}, failed (rho not found) {[w for w in results if results[w] is None]}')
    P(f'# components (lifts) with m|w|<=L: {sum(len(results[w]["rows"]) for w in okw)}; rows {len(out_rows)}; '
      f'hyperbolic clean rows {len(clean)}; distinct volumes {distinct}; '
      f'non-hyperbolic rows {sum(1 for r in out_rows if not r["hyperbolic"])}')
    P(f'# equal-volume pairs between rows: {cross}')
    P(f'# suspects/excluded: {suspects}')
    with open(f'summary_{surface}_L{L:g}.txt', 'w') as f:
        f.write('\n'.join(log))


if __name__ == '__main__':
    run(sys.argv[1], float(sys.argv[2]))
