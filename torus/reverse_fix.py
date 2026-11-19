"""Hyperbolic structures for lifts that fill S but on which SnapPy found none, from the reversed curve.

usage: python -I reverse_fix.py <surface> <workdir> <L> [procs]      (after harden.py, before merge_parts.py)

A closed geodesic gamma and its reverse gamma^-1 have isometric complements M_gamma and M_gamma^-1, by the
involution v -> -v of UT(S).  The reverse of the lift of w is a lift of the code word swap(reverse(w)), with the
same unoriented surface word and degree m.  For every lift that fills (exact test, filling_exact.py) but has
no positively oriented solution, the reversed code words are run in <workdir>_reversed (run_surface.py, then
harden.py), and where the reversed lift is hyperbolic its triangulation and volume are copied into the original
row, with status 'hyperbolic via the reversed curve (v -> -v)'.
"""
import sys, os, json, glob, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
surface, workdir, L = sys.argv[1], sys.argv[2], sys.argv[3]
procs = sys.argv[4] if len(sys.argv) > 4 else '2'
saved = sys.argv
sys.argv = [saved[0], surface, os.devnull]          # filling_exact reads its own argv
import filling_exact as fe
sys.argv = saved
import psl2
import surface_words as sw
from lifts import _cycles
from pipeline import _basis

G = fe.G
B = _basis(G)
POS = 'all tetrahedra positively oriented'


def reverse(w):
    v = w[::-1].translate(str.maketrans('XY', 'YX'))
    return min(v[i:] + v[:i] for i in range(len(v)))


def curves(w):
    """{(unoriented surface word, m): fills} for the lifts of w"""
    Mw = psl2.xy_matrix(w)
    out = {}
    for cy in _cycles(G.perm(Mw)):
        m = len(cy)
        P = psl2.I2
        for _ in range(m):
            P = psl2.mul(P, Mw)
        core = psl2.mul(psl2.mul(G.reps[cy[0]], P), psl2.inv(G.reps[cy[0]]))
        out[(sw.letters(sw.unoriented_canon(B.word(core))), m)] = fe.fills(w, cy[0], m)
    return out


todo = []                                           # (json path, row index, label, m, reversed word)
for f in sorted(glob.glob(os.path.join(workdir, '*.json'))):
    d = json.load(open(f))
    if not d.get('ok'):
        continue
    fills = None
    for i, r in enumerate(d['rows']):
        if r['hyperbolic']:
            continue
        fills = fills or curves(d['word'])
        lab = r.get('surface_word')
        ok = fills.get((lab, int(r['m']))) if lab else any(v for (l, m), v in fills.items() if m == int(r['m']))
        if ok:
            todo.append((f, i, lab, int(r['m']), reverse(d['word'])))
print(f'{len(todo)} lifts fill but have no hyperbolic structure', flush=True)
if not todo:
    sys.exit()

rdir = workdir.rstrip('/\\') + '_reversed'
wf = rdir + '_words.txt'
open(wf, 'w').write('\n'.join(sorted({t[4] for t in todo})) + '\n')
log = open(rdir + '.log', 'a')
subprocess.run([sys.executable, '-I', '-X', 'faulthandler', os.path.join(HERE, 'run_surface.py'), '--surface', surface,
                '--L', L, '--word-file', wf, '--workdir', rdir, '--out', rdir + '.csv', '--procs', procs,
                '--timeout', '3600'], stdout=log, stderr=subprocess.STDOUT, check=False)
subprocess.run([sys.executable, '-I', os.path.join(HERE, 'harden.py'), surface, rdir, procs], stdout=log,
               stderr=subprocess.STDOUT, check=False)

fixed = 0
for f, i, lab, m, rw in todo:
    rf = os.path.join(rdir, rw + '.json')
    if not os.path.exists(rf):
        print('  no result for the reverse of', os.path.basename(f)[:-5]); continue
    R = json.load(open(rf))
    match = [r for r in R.get('rows', []) if r['hyperbolic'] and int(r['m']) == m and (lab is None or r.get('surface_word') == lab)]
    if not match:
        print('  reverse of', os.path.basename(f)[:-5], 'is not hyperbolic either'); continue
    d = json.load(open(f))
    row, src = d['rows'][i], match[0]
    for k in ('volume', '_vol', '_tri', '_Kpos', 'solution_type'):
        row[k] = src[k]
    row['hyperbolic'] = True
    row['status'] = (row['status'] + '; ' if row['status'] != 'ok' else '') + f'hyperbolic via the reversed curve (v -> -v), {rw}'
    json.dump(d, open(f, 'w'))
    fixed += 1
    print(f'  {os.path.basename(f)[:-5]} {lab} m={m}: volume {src["volume"]} from {rw}')
print(f'fixed {fixed} of {len(todo)}')
