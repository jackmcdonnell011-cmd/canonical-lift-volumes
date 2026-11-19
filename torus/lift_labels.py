"""Exact labelling of the K-lifts in the cover: which curve on S each cusp of the cover over K is.

The cover of M = S^3 - (trefoil u K) is built by SnapPy from a permutation representation, with sheet x of
base tetrahedron t as cover tetrahedron t*d + x.  SnapPy glues sheets across a face by the permutation of
that face's generator in the unsimplified presentation (0 on faces of its spanning tree).  So for every
corner (t, v) of the cusp of K the sheets at that corner split into cover cusps exactly as the orbits of
the loops in the cusp torus based at (t, v), read as words in those face generators.  These loops are
computed here directly, so no corner has to be guessed.

Checks (else None is returned): the predicted split must equal SnapPy's at every corner of the cusp of K;
every convention that passes must give the same labels; and the labelled cusps with their degrees must be
exactly the curves of the coset cycles of M_w with their lengths, as a multiset.
"""
import psl2
from lifts import _cycles
import surface_words as sw


def _orbits(perms, d):
    parent = list(range(d))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for p in perms:
        for x in range(d):
            a, b = find(x), find(p[x])
            if a != b:
                parent[a] = b
    out = {}
    for x in range(d):
        out.setdefault(find(x), set()).add(x)
    return {frozenset(s) for s in out.values()}


def exact_labels(ml, G, B, C=None):
    """{cover cusp over K: unoriented cyclic surface word}, or None if the self-check fails."""
    M, gm, Gp = ml.M, ml.gm, ml.G
    d = G.d
    G0 = M.fundamental_group(simplify_presentation=False)
    om = [psl2.word_matrix(x, gm) for x in Gp.original_generators()]
    assert len(om) == len(G0.generators())
    M._choose_generators(False, False)
    info = M._choose_generators_info()
    bidx = M._get_cusp_indices_and_peripheral_curve_data()[0]
    if C is None:
        C = ml.cover(G).C
    cidx = C._get_cusp_indices_and_peripheral_curve_data()[0]
    K = ml.cuspK
    corners = [(t, v) for t in range(len(bidx)) for v in range(4) if bidx[t][v] == K]

    def edges(t, v):
        for f in range(4):
            if f == v:
                continue
            n, gl, g = info[t]['neighbors'][f], info[t]['gluings'][f], info[t]['generators'][f]
            yield (n, gl[v]), g

    def gmat(g, flip):
        if g == 0:
            return psl2.I2
        A = om[abs(g) - 1]
        return A if (g > 0) != flip else psl2.inv(A)

    def loops(root, flip):
        """matrices of loops in the cusp torus based at root (one per non-tree edge of a BFS tree)"""
        path = {root: psl2.I2}
        order, out = [root], []
        i = 0
        while i < len(order):
            a = order[i]; i += 1
            for b, g in edges(*a):
                e = psl2.mul(path[a], gmat(g, flip))
                if b not in path:
                    path[b] = e
                    order.append(b)
                else:
                    out.append(psl2.mul(e, psl2.inv(path[b])))
        return out

    def actual(t, v):
        part = {}
        for x in range(d):
            part.setdefault(cidx[t * d + x][v], set()).add(x)
        return part

    def label(L):
        """{cusp: unoriented surface word} from the loops L at every corner, or None if corners disagree"""
        labels = {}
        for (t, v) in corners:
            hyp = [A for A in L[(t, v)] if abs(psl2.tr(A)) > 2]
            if not hyp:
                continue
            g = min(hyp, key=lambda A: abs(psl2.tr(A)))
            pg = G.perm(g)
            for c, xs in actual(t, v).items():
                x = min(xs)
                p, y = 1, pg[x]
                while y != x:
                    y = pg[y]; p += 1
                gp = g
                for _ in range(p - 1):
                    gp = psl2.mul(gp, g)
                core = psl2.mul(psl2.mul(G.reps[x], gp), psl2.inv(G.reps[x]))
                if G.act(0, core) != 0:
                    return None
                u = sw.cyclic_canon(B.word(core))
                n = len(u)
                for q in range(1, n + 1):                 # primitive root of the cyclic word
                    if n % q == 0 and list(u[:q]) * (n // q) == list(u):
                        u = list(u[:q]); break
                lab = sw.letters(sw.unoriented_canon(u))
                if labels.setdefault(c, lab) != lab:
                    return None                           # two corners disagree: should never happen
        return labels

    # every convention (direction of the face generators, perm or inverse perm) that matches at every
    # corner is labelled; they must all give the same labels
    found = []
    for flip in (False, True):
        for inv in (False, True):
            L = {}
            for (t, v) in corners:
                mats = loops((t, v), flip)
                if inv:
                    mats = [psl2.inv(A) for A in mats]
                L[(t, v)] = mats
                pred = _orbits([G.perm(A) for A in mats], d)
                if pred != {frozenset(s) for s in actual(t, v).values()}:
                    break
            else:
                found.append(label(L))
    if not found or any(lab is None or lab != found[0] for lab in found):
        return None
    labels = found[0]

    # the labelled cusps, with their degrees (sheets per corner), must be exactly the curves of the coset
    # cycles of M_w with their lengths, as a multiset: on a reciprocal w a curve and its reverse both occur
    from collections import Counter
    Mw = psl2.xy_matrix(ml.w)
    want = Counter()
    for cy in _cycles(G.perm(Mw)):
        P = psl2.I2
        for _ in range(len(cy)):
            P = psl2.mul(P, Mw)
        core = psl2.mul(psl2.mul(G.reps[cy[0]], P), psl2.inv(G.reps[cy[0]]))
        want[(sw.letters(sw.unoriented_canon(B.word(core))), len(cy))] += 1
    t0, v0 = corners[0]
    deg = {c: len(xs) for c, xs in actual(t0, v0).items()}
    if Counter((labels.get(c), m) for c, m in deg.items()) != want:
        return None
    return labels
