"""Exact filling test for a lift of a modular geodesic to S = H/Gamma, independent of SnapPy.

usage: python -I filling_exact.py <surface> <results.csv> [--only-open]

S is cut into the ideal Farey triangles.  The lift crosses one Farey edge per letter (see
filling_check.py).  In the frame where that edge is (0, oo), the crossing point is i*t with
t^2 = b/c for the matrix [[a,b],[c,d]] of the rotated word (the axis has endpoints r1 < 0 < r2 and
t^2 = -r1*r2 = b/c), so all positions are exact rationals.
Inside a triangle the curve is a set of geodesic chords; two boundary pieces lie in the same
complementary face iff no chord separates them (half-planes are convex).  Gluing the faces across
the pieces of Farey edges gives the components of (closed-up S) minus the curve, and
chi = faces - edge pieces + cusps.  The lift fills S iff every component has chi = 1 and contains
at most one cusp (a disk or a once-punctured disk).  Faces touching no Farey edge are disks.
"""
import csv, sys, os
csv.field_size_limit(2**31 - 1)      # status fields can list many equal-volume partners
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import psl2
from gn_volumes import SURFACES
from lifts import _cycles

R = psl2.mul(psl2.S, psl2.X)
G = None


def setup(Gamma):
    """Prepare fills() for the surface H/Gamma (Gamma: a subgroups.CosetAction or a surface name)."""
    global G, d, pS, pR, pX, pY, pT, cusp_of, tri_of, side_k
    G = SURFACES[Gamma]() if isinstance(Gamma, str) else Gamma
    d = G.d
    pS, pR, pX, pY = (G.perm(M) for M in (psl2.S, R, psl2.X, psl2.Y))
    pT = G.perm(psl2.X)
    cusp_of = {}
    for k, cy in enumerate(_cycles(pT)):
        for c in cy:
            cusp_of[c] = k
    # triangle frames: ccw boundary order of the triangle of dart d0 is side d0, d0.R^2, d0.R, each with t increasing
    tri_of, side_k = {}, {}
    for d0 in range(d):
        if d0 in tri_of:
            continue
        for k, x in enumerate((d0, pR[pR[d0]], pR[d0])):
            tri_of[x], side_k[x] = d0, k


# callers may set sys.argv = [prog, surface, path] before importing
if __name__ != '__main__' and len(sys.argv) >= 2 and sys.argv[1] in SURFACES:
    setup(sys.argv[1])


def axis_key(A):
    """Order of the crossing of axis(A) with (0, oo): height t^2 = b/c, then minus the centre
    (a-d)/(2c).  The second entry breaks ties (curve crossing itself on the edge) by pushing the
    crossing slightly into the triangle on the Re > 0 side; both triangles see the same order."""
    a, b, c, dd = A
    return (Fraction(b, c), -Fraction(a - dd, 2 * c))


def fills(w, c0, m):
    n = len(w)
    rotM = [psl2.xy_matrix(w[j:] + w[:j]) for j in range(n)]
    Sinv = psl2.inv(psl2.S)
    darts = []
    c = c0
    for s in range(m * n):
        darts.append(c)
        c = (pX if w[s % n] == 'X' else pY)[c]
    assert c == c0
    N = len(darts)
    # crossing s lies on edge e_s = min(c_s, c_s.S); order the crossings along each edge in its frame
    on_edge = {}
    for s in range(N):
        x = darts[s]
        e = min(x, pS[x])
        A = rotM[s % n] if x == e else psl2.mul(psl2.mul(psl2.S, rotM[s % n]), Sinv)
        on_edge.setdefault(e, []).append((axis_key(A), s))
    rank = {}
    for e, lst in on_edge.items():
        lst.sort()
        if len(set(k for k, _ in lst)) != len(lst):
            return None     # two crossings with the same axis: should not happen
        for i, (_, s) in enumerate(lst):
            rank[s] = i
    def place(x, i, k):
        """position along side x (t increasing in x's frame) of rank i on an edge with k points"""
        return i if x == min(x, pS[x]) else k - 1 - i
    def key(x, p):          # position around the triangle boundary
        return (side_k[x], p)
    chords = {}
    for s in range(N):
        s1 = (s + 1) % N
        a_dart, b_dart = pS[darts[s]], darts[s1]
        assert tri_of[a_dart] == tri_of[b_dart]
        ka = key(a_dart, place(a_dart, rank[s], len(on_edge[min(a_dart, pS[a_dart])])))
        kb = key(b_dart, place(b_dart, rank[s1], len(on_edge[min(b_dart, pS[b_dart])])))
        chords.setdefault(tri_of[a_dart], []).append((min(ka, kb), max(ka, kb)))
    # edge pieces: edge e (frame e) is cut at its points into pieces 0..k; piece i, seen from dart x
    parent = {}
    def find(u):
        parent.setdefault(u, u)
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u
    def union(u, v):
        parent[find(u)] = find(v)
    faces, pieces = set(), 0
    piece_face = {}
    for e in range(d):
        if e != min(e, pS[e]):
            continue
        k = len(on_edge.get(e, []))
        for i in range(k + 1):
            # piece i lies between points i-1 and i in e's frame
            pieces += 1
            fs = []
            for x, p in ((e, Fraction(2 * i - 1, 2)), (pS[e], Fraction(2 * k - 2 * i - 1, 2))):
                T = tri_of[x]
                kq = key(x, p)
                sig = tuple(lo < kq < hi for lo, hi in chords.get(T, []))
                f = (T, sig)
                faces.add(f)
                fs.append(f)
            union(fs[0], fs[1])
            piece_face[(e, i)] = fs
    # cusps: the corner at t = oo of side x is cusp_of[x]; its face is that of the last piece of side x
    cusp_face = {}
    for x in range(d):
        e = min(x, pS[x])
        i = len(on_edge.get(e, [])) if x == e else 0
        f = piece_face[(e, i)][0 if x == e else 1]
        cusp_face.setdefault(cusp_of[x], []).append(f)
    for k, fs in cusp_face.items():
        for f in fs[1:]:
            union(fs[0], f)
    chi, ncusp = {}, {}
    for f in faces:
        r = find(f)
        chi[r] = chi.get(r, 0) + 1
    for (e, i), fs in piece_face.items():
        r = find(fs[0])
        chi[r] -= 1
    for k, fs in cusp_face.items():
        r = find(fs[0])
        chi[r] += 1
        ncusp[r] = ncusp.get(r, 0) + 1
    return all(chi[r] == 1 and ncusp.get(r, 0) <= 1 for r in chi)


if __name__ == '__main__':
    surface, path = sys.argv[1], sys.argv[2]
    setup(surface)
    cache = {}
    def lifts(w, m):
        if w not in cache:
            cache[w] = [(len(cy), fills(w, cy[0], len(cy))) for cy in _cycles(G.perm(psl2.xy_matrix(w)))]
        return [f for mm, f in cache[w] if mm == m]


    rows = list(csv.DictReader(open(path, encoding='utf-8')))
    groups = {}
    for r in rows:
        groups.setdefault((r['word'], int(r['m'])), []).append(r)
    agree = disagree = undecided = 0
    bad = []
    newfill = []
    for (w, m), rs in groups.items():
        fl = lifts(w, m)
        if len(rs) == 1 and len(fl) > 1 and len(set(fl)) == 1:
            fl = fl[:1]
        if None in fl or len(rs) != len(fl):
            undecided += len(rs)
            continue
        nhyp = sum(str(r['hyperbolic']) == 'True' for r in rs)
        nfill = sum(fl)
        if nhyp == nfill:
            agree += len(rs)
        else:
            disagree += len(rs)
            bad.append((w, m, f'hyperbolic {nhyp}/{len(rs)}', f'filling {nfill}/{len(fl)}'))
    print(f'{surface} {os.path.basename(path)}: {len(rows)} rows')
    print(f'  SnapPy hyperbolic == exact filling: {agree}   disagree: {disagree}   undecided: {undecided}')
    for b in bad[:20]:
        print('   ', b)
    # The counts above are per (word, m).  Rows with a surface-word label (exact labels) are also
    # compared one by one: the label picks out the coset cycle, whose filling is tested directly.
    import surface_words as sw
    B = sw.Basis(G)
    curve = {}
    def cycle_fills(w):
        if w not in curve:
            Mw, out = psl2.xy_matrix(w), {}
            for cy in _cycles(G.perm(Mw)):
                P = psl2.I2
                for _ in range(len(cy)):
                    P = psl2.mul(P, Mw)
                core = psl2.mul(psl2.mul(G.reps[cy[0]], P), psl2.inv(G.reps[cy[0]]))
                out.setdefault((sw.letters(sw.unoriented_canon(B.word(core))), len(cy)), set()).add(fills(w, cy[0], len(cy)))
            curve[w] = out
        return curve[w]
    same = diff = nolabel = 0
    for r in rows:
        f = cycle_fills(r['word']).get((r.get('surface_word') or '', int(r['m'])))
        if not f or len(f) != 1:
            nolabel += 1
            continue
        if (str(r['hyperbolic']) == 'True') == next(iter(f)):
            same += 1
        else:
            diff += 1
            print('    per-lift disagreement:', r['word'], r['m'], r.get('surface_word'), r['hyperbolic'])
    print(f'  per lift (labelled rows): agree {same}   disagree {diff}   no usable label {nolabel}')
