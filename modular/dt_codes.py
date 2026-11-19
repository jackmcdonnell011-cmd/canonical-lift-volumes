"""Code words of closed geodesics on the modular surface, and their DT codes on the Lorenz template.

A closed geodesic on Mod = PSL2(Z)\\H^2 is recorded by its cyclic code word w in the letters X (left turn)
and Y (right turn).  Every geodesic is listed once by taking only primitive words (not a power of a shorter
word) in their lexicographically smallest rotation, i.e. the binary Lyndon words.  The words X and Y are
parabolic and are left out.

get_dt converts a code word into the Dowker-Thistlethwaite code of the two-component link formed by the
Lorenz knot of w and the trefoil, using Williams' description of the Lorenz template.  It is the function
get_dt of notebooks/Modular_DT_Codes.ipynb, unchanged apart from comments.

usage:  python dt_codes.py <max_length> <out.csv>      writes word,dt_code for all words of length 2..max_length
"""
import csv
import sys
from collections import defaultdict


def is_primitive(w):
    n = len(w)
    return not any(w[:d] * (n // d) == w for d in range(1, n // 2 + 1) if n % d == 0)


def canonical(w):
    """the lexicographically smallest rotation of w (X < Y)"""
    n = len(w)
    return min((w + w)[i:i + n] for i in range(n))


def code_words(max_length, min_length=2):
    """All primitive canonical code words of length min_length..max_length (Duval's algorithm for Lyndon
    words), in order of length and then lexicographically."""
    by_length = defaultdict(list)
    w = [-1]
    while w:
        w[-1] += 1
        if len(w) >= min_length:
            by_length[len(w)].append(''.join('XY'[c] for c in w))
        m = len(w)
        while len(w) < max_length:
            w.append(w[len(w) - m])
        while w and w[-1] == 1:
            w.pop()
    for n in range(min_length, max_length + 1):
        yield from sorted(by_length[n])


def get_dt(w):
    """DT code of (Lorenz knot of w) + trefoil, as the string '[(...),(...)]'.  None if w is not primitive
    and canonical."""
    n = len(w)
    # filters: primitive and canonical rotation
    if any(w[:d] * (n // d) == w for d in range(1, n // 2 + 1) if n % d == 0):
        return None
    if w != min((w + w)[i:i + n] for i in range(n)):
        return None

    # position of each strand on the branch line = rank of the corresponding cyclic shift
    wt = [0] * n
    for r, (_, i) in enumerate(sorted([(w[i:] + w[:i], i) for i in range(n)])):
        wt[i] = r
    strands = [{'s': wt[i], 'e': wt[(i + 1) % n]} for i in range(n)]

    # crossings between strands, ordered along each strand
    cmap, cid = defaultdict(list), 1
    for i in range(n):
        for j in range(i + 1, n):
            s1, s2 = strands[i], strands[j]
            d1, d2 = s1['e'] - s1['s'], s2['e'] - s2['s']
            if d1 == d2:
                continue
            t = (s2['s'] - s1['s']) / (d1 - d2)
            if 0 < t < 1:
                c = {'id': cid, 't': t}
                cmap[i].append(c)
                cmap[j].append(c)
                cid += 1
    for i in range(n):
        cmap[i].sort(key=lambda x: x['t'])

    # traverse the knot, labelling its crossings and its passes over/under the trefoil
    top, bot, xl, ctr = defaultdict(list), defaultdict(list), defaultdict(list), 3
    curr = wt.index(0)
    top[wt[curr]].append(1)
    bot[wt[curr]].append(-2)
    for _ in range(n):
        s = strands[curr]
        sign = 1 if s['e'] > s['s'] else -1
        for c in cmap[curr]:
            xl[c['id']].append(ctr * sign)
            ctr += 1
        ew = s['e']
        if ew == 0:
            break
        top[ew].append(ctr)
        ctr += 1
        bot[ew].append(-ctr)
        ctr += 1
        curr = wt.index(ew)

    # traverse the trefoil
    mv, mx = wt[(wt.index(max(wt)) - 1 + n) % n], max(wt)

    def fill(r, d, s):
        nonlocal ctr
        for v in r:
            if len(d[v]) < 2:
                d[v].append(-ctr if s else ctr)
                ctr += 1
    fill(range(mv + 1), bot, 0)
    t1 = list(range(ctr, ctr + 3))
    ctr += 3
    fill(range(mv + 1, mx + 1), top, 1)
    fill(range(mx, mv, -1), bot, 0)
    t2 = list(range(ctr, ctr + 3))
    ctr += 3
    fill(range(mv, -1, -1), top, 1)

    # pair the labels at each crossing and write out the DT code
    pairs = [top[x] for x in wt] + [bot[x] for x in wt] + list(xl.values())
    tf = [t1[0], -t1[1], t1[2], -t2[0], t2[1], -t2[2]]
    pairs += [[tf[i], tf[i + 3]] for i in range(3)]
    vp = [p for p in pairs if len(p) == 2]
    if not vp:
        return "([],[])"
    dt = [0] * ((max(max(abs(x) for x in p) for p in vp) // 2) + 2)
    for o, e in vp:
        if abs(o) % 2 == 0:
            o, e = e, o
        if (idx := (abs(o) - 1) // 2) < len(dt):
            dt[idx] = e
    while dt and dt[-1] == 0:
        dt.pop()
    s = dt.index(-2)
    return f"[{tuple(dt[:s])},{tuple(dt[s:])}]"


if __name__ == '__main__':
    max_length, out = int(sys.argv[1]), sys.argv[2]
    with open(out, 'w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['word', 'dt_code'])
        for w in code_words(max_length):
            wr.writerow([w, get_dt(w)])
