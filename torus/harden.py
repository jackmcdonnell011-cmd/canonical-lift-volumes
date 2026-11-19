"""Retry SnapPy's hyperbolic structure where the exact filling test says it must exist.

usage: python -I harden.py <surface> <workdir> [procs]

For every (word, m) where fewer lifts came out hyperbolic than the exact test (filling_exact.py) says
fill, every non-hyperbolic lift of that (word, m) is retried: up to 300 random retriangulations, with
and without simplifying first.  A lift that reaches 'all tetrahedra positively oriented' is updated in
its word JSON (volume, solution type, triangulation) and its status gets 'hyperbolic after retry'.
Rerun merge_parts.py afterwards.  Prints any (word, m) that still does not match.
"""
import sys, os, json, glob
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
POS = 'all tetrahedra positively oriented'


def retry(item):
    import snappy
    key, idx, tri = item
    for simplify in (False, True):
        F = snappy.Manifold(tri)
        if simplify:
            F.simplify()
        for _ in range(150):
            if F.solution_type() == POS:
                return key, idx, float(F.volume()), F._to_string()
            F.randomize()
    return key, idx, None, None


def main():
    surface, workdir = sys.argv[1], sys.argv[2]
    procs = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    saved = sys.argv
    sys.argv = [saved[0], surface, os.devnull]
    import filling_exact as fe
    sys.argv = saved
    import psl2
    from lifts import _cycles
    todo, files = [], {}
    for f in sorted(glob.glob(os.path.join(workdir, '*.json'))):
        d = json.load(open(f))
        if not d.get('ok'):
            continue
        w = d['word']
        groups = {}
        for i, r in enumerate(d['rows']):
            groups.setdefault(int(r['m']), []).append(i)
        nfill = {}
        for cy in _cycles(fe.G.perm(psl2.xy_matrix(w))):
            nfill.setdefault(len(cy), []).append(fe.fills(w, cy[0], len(cy)))
        for m, idxs in groups.items():
            fl = nfill.get(m, [])
            if len(idxs) == 1 and len(fl) > 1 and len(set(fl)) == 1:
                fl = fl[:1]          # normal Gamma: one row stands for all lifts
            nhyp = sum(str(d['rows'][i]['hyperbolic']) == 'True' for i in idxs)
            if nhyp < sum(fl):
                files[f] = d
                for i in idxs:
                    if str(d['rows'][i]['hyperbolic']) != 'True':
                        todo.append(((f, m), i, d['rows'][i]['_tri']))
    print(f'{len(todo)} lifts to retry in {len(files)} words', flush=True)
    fixed = 0
    with Pool(procs) as P:
        for (f, m), i, v, tri in P.imap_unordered(retry, todo):
            if v is None:
                continue
            r = files[f]['rows'][i]
            r.update(volume=round(v, 10), _vol=v, _tri=tri, hyperbolic=True, solution_type=POS,
                     status=(r['status'] + '; ' if r['status'] != 'ok' else '') + 'hyperbolic after retry')
            fixed += 1
            print(f"{files[f]['word']} m={m} cusp {r['lift_cusp']}: vol {v:.6f}", flush=True)
    for f, d in files.items():
        json.dump(d, open(f + '.tmp', 'w'))
        os.replace(f + '.tmp', f)
    print(f'fixed {fixed} of {len(todo)}', flush=True)


if __name__ == '__main__':
    main()
