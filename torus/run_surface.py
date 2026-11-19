"""Single entry point: canonical-lift volumes on one surface, one chunk of words.

  python -I run_surface.py --surface torus11 --L 22 --word-file words.txt --workdir work --out results.csv --procs 4

* Words: every primitive modular code word w (letters X,Y, both present, canonical
  rotation) with |w| <= L that has at least one lift with m*|w| <= L.  Sorted
  deterministically; chunk i/k takes words i-1, i-1+k, i-1+2k, ...
* Every word is processed with all checks (see pipeline.process_word) and saved to
  <workdir>/<word>.json as soon as it finishes, so a rerun resumes where it stopped.
* rho (the Ghys map) is read from rho_cache/ if present, else solved and cached.
* Writes --out (CSV) and --out with .summary.txt.
"""
import os
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import argparse
import csv
import json
import time
import traceback

POS = 'all tetrahedra positively oriented'
_STATE = {}


def _setup(surface, L, cache):
    import psl2
    from gn_volumes import SURFACES
    from pipeline import trefoil_cover, cover_fingerprint, normaliser_autos
    G = SURFACES[surface]()
    sig = G.signature()
    autos = normaliser_autos(G)
    _STATE.update(G=G, sig=sig, normal=(len(autos) == G.d), Dfp=cover_fingerprint(trefoil_cover(G)),
                  L=L, cache=cache, surface=surface)


def word_list(surface, L):
    import psl2
    from gn_volumes import SURFACES
    from lifts import _cycles
    from build_rho_cache import words_upto
    G = SURFACES[surface]()
    out = []
    for w in words_upto(L):
        cyc = _cycles(G.perm(psl2.xy_matrix(w)))
        if min(len(c) for c in cyc) * len(w) <= L:
            out.append(w)
    return out


def filling_words(surface, L, words):
    """Words with at least one lift of length m|w| <= L that fills S.  Non-filling lifts only ever
    give non-hyperbolic complements, so for volumes nothing is lost."""
    import psl2
    from lifts import _cycles
    saved = sys.argv
    sys.argv = [saved[0], surface, os.devnull]     # filling_exact reads its own argv
    try:
        import filling_exact as fe
    finally:
        sys.argv = saved
    out = []
    for w in words:
        for cy in _cycles(fe.G.perm(psl2.xy_matrix(w))):
            if len(cy) * len(w) <= L and fe.fills(w, cy[0], len(cy)):
                out.append(w)
                break
    return out


def process(w):
    """Run one word; return a JSON-able dict (never raises)."""
    from pipeline import process_word, cusp_iso
    st = _STATE
    G, L = st['G'], st['L']
    t0 = time.time()
    if st.get('workdir'):
        # start marker, read by the watchdog in run_pool
        open(os.path.join(st['workdir'], w + '.running'), 'w').write(str(t0))
    try:
        r = process_word(w, G, L, st['Dfp'], st['sig'], cache_dir=st['cache'])
    except Exception as e:
        return dict(word=w, ok=False, error=f'{type(e).__name__}: {e}', trace=traceback.format_exc()[-800:],
                    seconds=round(time.time() - t0, 1))
    rows = r['rows']
    checks = {k: bool(v) for k, v in r['checks'].items()}
    allok = all(checks.values())
    groups = []
    if rows:
        if st['normal']:
            status = []
            m = rows[0]['m']
            if len(set(r['degrees'])) != 1 or len(r['degrees']) != G.d // m or len(rows) != G.d // m:
                status.append(f'SUSPECT: degrees {r["degrees"]}, {len(rows)} lifts, expected index/m equal lifts')
            hypw = all(lf['solution'] == POS for lf in rows)
            if hypw:
                vols = [lf['volume'] for lf in rows]
                if max(vols) - min(vols) > 1e-8:
                    status.append('SUSPECT: lifts of one word have different volumes')
                for lf in rows[1:]:
                    if cusp_iso(rows[0]['F'], rows[0]['Kpos'], lf['F'], lf['Kpos']) is not True:
                        status.append(f'SUSPECT: lift cusp {lf["cusp"]} not cusp-isometric to cusp {rows[0]["cusp"]}')
            groups.append((rows[0], len(rows), status, hypw))
        else:
            for lf in rows:
                groups.append((lf, 1, [], lf['solution'] == POS))
    out_rows = []
    ck = ' '.join(f'{k}={"pass" if v else "FAIL"}' for k, v in checks.items())
    for lf, mult, status, hypc in groups:
        if not allok:
            status = status + ['EXCLUDED: a check failed']
        out_rows.append(dict(
            word=w, modular_word_length=len(w), m=lf['m'], surface_word_length=lf['m'] * len(w),
            hyperbolic_length_on_S=round(lf['length'], 10),
            volume=round(lf['volume'], 10) if hypc else '',
            multiplicity=mult, gamma_trace=lf['trace'], checks=ck,
            gamma_element=str(lf['gamma']) + ('' if (lf['gamma_pinned'] or st['normal']) else ' (one of the equal-m components)'),
            hyperbolic=hypc, solution_type=lf['solution'],
            status='; '.join(status) if status else 'ok',
            lift_cusp=lf['cusp'], surface_word=lf.get('surface_word'), _tri=lf['F']._to_string(), _Kpos=lf['Kpos'], _vol=lf['volume']))
    return dict(word=w, ok=True, checks=checks, degrees=r['degrees'], g_cover=r['g_cover'], n_cover=r['n_cover'],
                solver=r['solver'], label_method=r.get('label_method'), rows=out_rows, seconds=round(time.time() - t0, 1))


def _worker_init(surface, L, cache, workdir=None):
    _setup(surface, L, cache)
    _STATE['workdir'] = workdir


def save(workdir, res):
    fn = os.path.join(workdir, res['word'] + '.json')
    json.dump(res, open(fn + '.tmp', 'w'))
    os.replace(fn + '.tmp', fn)
    if res['ok']:
        print(f"{res['word']}: checks {'PASS' if all(res['checks'].values()) else 'FAIL ' + str(res['checks'])} "
              f"rows {len(res['rows'])} [{res['seconds']}s]", flush=True)
    else:
        print(f"{res['word']}: SKIPPED ({res['error']}) [{res['seconds']}s]", flush=True)


def run_pool(a, todo, workdir):
    """Pool run with a per-word time limit.  A word running longer than a.timeout seconds is saved
    as skipped (ok=False), the pool is restarted and the unfinished words are run again."""
    from multiprocessing import Pool
    remaining = list(todo)
    while remaining:
        for f in os.listdir(workdir):
            if f.endswith('.running'):
                os.remove(os.path.join(workdir, f))
        pool = Pool(a.procs, initializer=_worker_init, initargs=(a.surface, a.L, a.cache, workdir))
        pending = {w: pool.apply_async(process, (w,)) for w in remaining}
        stuck = None
        while pending and stuck is None:
            time.sleep(2)
            for w, ar in list(pending.items()):
                if ar.ready():
                    save(workdir, ar.get())
                    del pending[w]
            for w in pending:
                mk = os.path.join(workdir, w + '.running')
                if os.path.exists(mk) and time.time() - float(open(mk).read() or time.time()) > a.timeout:
                    stuck = w
                    break
        pool.terminate()
        pool.join()
        if stuck:
            save(workdir, dict(word=stuck, ok=False, error=f'timeout: no result after {a.timeout}s',
                               seconds=a.timeout))
            del pending[stuck]
        remaining = list(pending)
    for f in os.listdir(workdir):
        if f.endswith('.running'):
            os.remove(os.path.join(workdir, f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--surface', required=True)
    ap.add_argument('--L', type=int, required=True, help='max surface word length m*|w|')
    ap.add_argument('--chunk', default='1/1', help='i/k, 1-based')
    ap.add_argument('--out', required=True)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--workdir', default=None)
    ap.add_argument('--cache', default=os.path.join(HERE, 'rho_cache'))
    ap.add_argument('--timeout', type=float, default=2400, help='seconds per word before it is skipped (procs > 1)')
    ap.add_argument('--filling-only', action='store_true',
                    help='only words with a lift (m|w| <= L) that fills S (exact test, filling_exact.py)')
    ap.add_argument('--word-file', default=None,
                    help='run exactly these modular words (one per line) instead of all words up to L')
    a = ap.parse_args()
    i, k = map(int, a.chunk.split('/'))
    workdir = a.workdir or os.path.join(HERE, f'work_{a.surface}_L{a.L}')
    os.makedirs(workdir, exist_ok=True)
    if a.word_file:
        words = sorted(set(open(a.word_file).read().split()))
    else:
        words = word_list(a.surface, a.L)
    if a.filling_only:
        words = filling_words(a.surface, a.L, words)
    words = words[i - 1::k]
    todo = [w for w in words if not os.path.exists(os.path.join(workdir, w + '.json'))]
    print(f'surface {a.surface} L={a.L} chunk {i}/{k}: {len(words)} words, {len(todo)} to do', flush=True)
    if todo:
        if a.procs > 1:
            run_pool(a, todo, workdir)
        else:
            _setup(a.surface, a.L, a.cache)
            for w in todo:
                save(workdir, process(w))
    # assemble this chunk
    rows, skipped, failed, suspect = [], [], [], []
    for w in words:
        res = json.load(open(os.path.join(workdir, w + '.json')))
        if not res['ok']:
            skipped.append(w)
            continue
        if not all(res['checks'].values()):
            failed.append(w)
        for r in res['rows']:
            rows.append({kk: vv for kk, vv in r.items() if not kk.startswith('_')})
            if 'SUSPECT' in r['status']:
                suspect.append(w)
    cols = ['word', 'modular_word_length', 'm', 'surface_word_length', 'hyperbolic_length_on_S', 'volume',
            'multiplicity', 'gamma_trace', 'checks', 'gamma_element', 'hyperbolic', 'solution_type', 'status', 'lift_cusp', 'surface_word']
    with open(a.out, 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=cols)
        wr.writeheader()
        wr.writerows(rows)
    summ = (f'surface={a.surface} L={a.L} chunk={a.chunk}\n'
            f'words={len(words)} rows={len(rows)} hyperbolic_rows={sum(1 for r in rows if r["hyperbolic"])}\n'
            f'skipped={len(skipped)} {skipped}\n'
            f'failed_checks={len(failed)} {failed}\n'
            f'suspect={len(set(suspect))} {sorted(set(suspect))}\n'
            f'workdir={workdir}\n')
    open(os.path.splitext(a.out)[0] + '.summary.txt', 'w').write(summ)
    print(summ, flush=True)


if __name__ == '__main__':
    main()
