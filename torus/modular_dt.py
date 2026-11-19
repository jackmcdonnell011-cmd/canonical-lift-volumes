# get_dt from notebooks/Modular_DT_Codes.ipynb (the same function as modular/dt_codes.py)
from collections import defaultdict

def get_dt(w):
    n = len(w)
    # Check Primitive & Canonical
    if any(w[:d]*(n//d)==w for d in range(1, n//2+1) if n%d==0): return None
    if w != min((w+w)[i:i+n] for i in range(n)): return None

    # Logic: Weights, Strands, Intersections
    wt = [0]*n
    for r, (_, i) in enumerate(sorted([(w[i:]+w[:i], i) for i in range(n)])): wt[i] = r
    strands = [{'s': wt[i], 'e': wt[(i+1)%n]} for i in range(n)]
    cmap, cid = defaultdict(list), 1
    for i in range(n):
        for j in range(i+1, n):
            s1, s2 = strands[i], strands[j]
            d1, d2 = s1['e']-s1['s'], s2['e']-s2['s']
            if d1 == d2: continue
            t = (s2['s']-s1['s']) / (d1-d2)
            if 0 < t < 1:
                c = {'id': cid, 't': t}; cmap[i].append(c); cmap[j].append(c); cid += 1
    for i in range(n): cmap[i].sort(key=lambda x: x['t'])

    # Logic: Passes & Labels
    top, bot, xl, ctr = defaultdict(list), defaultdict(list), defaultdict(list), 3
    curr = wt.index(0); top[wt[curr]].append(1); bot[wt[curr]].append(-2)
    for _ in range(n):
        s = strands[curr]; sign = 1 if s['e'] > s['s'] else -1
        for c in cmap[curr]: xl[c['id']].append(ctr * sign); ctr += 1
        ew = s['e']
        if ew == 0: break
        top[ew].append(ctr); ctr+=1; bot[ew].append(-ctr); ctr+=1; curr = wt.index(ew)

    # Logic: Padding & Trefoil
    mv, mx = wt[(wt.index(max(wt))-1+n)%n], max(wt)
    def fill(r, d, s):
        nonlocal ctr
        for v in r:
            if len(d[v])<2: d[v].append(-ctr if s else ctr); ctr+=1
    fill(range(mv+1), bot, 0); t1=list(range(ctr, ctr+3)); ctr+=3
    fill(range(mv+1, mx+1), top, 1); fill(range(mx, mv, -1), bot, 0)
    t2=list(range(ctr, ctr+3)); ctr+=3; fill(range(mv, -1, -1), top, 1)

    # Logic: Assembly
    pairs = [top[x] for x in wt] + [bot[x] for x in wt] + list(xl.values())
    tf = [t1[0], -t1[1], t1[2], -t2[0], t2[1], -t2[2]]
    pairs += [[tf[i], tf[i+3]] for i in range(3)]
    vp = [p for p in pairs if len(p)==2]
    if not vp: return "([],[])"

    dt = [0] * ((max(max(abs(x) for x in p) for p in vp)//2)+2)
    for o, e in vp:
        if abs(o)%2==0: o, e = e, o
        if (idx := (abs(o)-1)//2) < len(dt): dt[idx] = e
    while dt and dt[-1] == 0: dt.pop()
    s = dt.index(-2)
    return f"[{tuple(dt[:s])},{tuple(dt[s:])}]"

