"""Rows for words whose Ghys map rho was not found, from the reversed curve.

usage: python -I reverse_rho.py <surface> <workdir> <L> [procs]      (after the rho rerun, before harden.py)

A closed geodesic gamma and its reverse gamma^-1 have isometric complements M_gamma and M_gamma^-1, by the
involution v -> -v of UT(S) (see reverse_fix.py).  The reverse of a lift of w is a lift of the code word
swap(reverse(w)), with the same unoriented surface word and degree m.  For a word w skipped with 'rho not
found', each row of the reversed word rw (from <workdir>, or run in <workdir>_reversed if absent) gives the row
of w with the same (surface word, m); the Gamma element and trace are recomputed from w's own coset cycle, and
the status records the source.  The original failure is kept as <word>.json.norho.
"""
import sys, os, json, glob, shutil, subprocess

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


def reverse(w):
    v = w[::-1].translate(str.maketrans('XY', 'YX'))
    return min(v[i:] + v[:i] for i in range(len(v)))


def mpow(A, k):
    P = psl2.I2
    for _ in range(k):
        P = psl2.mul(P, A)
    return P


def gam_by_label(w):
    """{(unoriented surface word, m): Gamma element} for the lifts of w"""
    Mw = psl2.xy_matrix(w)
    out = {}
    for cy in _cycles(G.perm(Mw)):
        m = len(cy)
        gam = psl2.mul(psl2.mul(G.reps[cy[0]], mpow(Mw, m)), psl2.inv(G.reps[cy[0]]))
        out.setdefault((sw.letters(sw.unoriented_canon(B.word(gam))), m), gam)
    return out


todo = []
for f in sorted(glob.glob(os.path.join(workdir, '*.json'))):
    d = json.load(open(f))
    if not d.get('ok') and 'rho not found' in str(d.get('error', '')):
        todo.append((f, d['word'], reverse(d['word'])))
print(f'{len(todo)} words without rho', flush=True)
if not todo:
    sys.exit()

missing = [rw for f, w, rw in todo if not json.load(open(os.path.join(workdir, rw + '.json'))).get('ok')] \
    if all(os.path.exists(os.path.join(workdir, t[2] + '.json')) for t in todo) else [t[2] for t in todo]
rdir = workdir.rstrip('/\\') + '_reversed'
if missing:
    wf = rdir + '_rho_words.txt'
    open(wf, 'w').write('\n'.join(sorted(set(missing))) + '\n')
    log = open(rdir + '.log', 'a')
    subprocess.run([sys.executable, '-I', '-X', 'faulthandler', os.path.join(HERE, 'run_surface.py'), '--surface',
                    surface, '--L', L, '--word-file', wf, '--workdir', rdir, '--out', rdir + '_rho.csv', '--procs',
                    procs, '--timeout', '3600'], stdout=log, stderr=subprocess.STDOUT, check=False)

fixed = 0
for f, w, rw in todo:
    rf = os.path.join(workdir, rw + '.json')
    R = json.load(open(rf)) if os.path.exists(rf) else {}
    if not R.get('ok'):
        rf = os.path.join(rdir, rw + '.json')
        R = json.load(open(rf)) if os.path.exists(rf) else {}
    if not R.get('ok') or w == rw:
        print('  no usable result for the reverse of', w); continue
    gl = gam_by_label(w)
    rows = []
    for src in R['rows']:
        lab, m = src.get('surface_word'), int(src['m'])
        gam = gl.get((lab, m))
        if gam is None:
            break
        row = dict(src)
        row.update(word=w, gamma_element=str(tuple(gam)), gamma_trace=str(abs(psl2.tr(gam))),
                   checks='via the reversed curve', status=f'rho not found; row from the reversed curve (v -> -v), {rw}')
        rows.append(row)
    if len(rows) != len(R['rows']) or set(gl) != set(gam_by_label(rw)):     # same curves, reversed
        print('  labels of', w, 'and', rw, 'do not match'); continue
    shutil.move(f, f + '.norho')
    D = {k: v for k, v in R.items() if k not in ('rows', 'word', 'checks')}
    D.update(word=w, ok=True, checks={}, label_method='reversed curve', solver='none (reversed curve ' + rw + ')',
             rows=rows)
    json.dump(D, open(f, 'w'))
    fixed += 1
    print(f'  {w}: {len(rows)} rows from {rw}, volumes {sorted({r["volume"] for r in rows})}')
print(f'fixed {fixed} of {len(todo)}')
