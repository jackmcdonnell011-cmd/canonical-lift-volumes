"""Verified computations (interval arithmetic) for the once-punctured torus data; run inside SageMath with SnapPy.

usage:  sage -python verify_sage.py <workdir> <out.csv> [procs]      (workdir = work_N<N> from run_torus.sh)

Each hyperbolic manifold is re-solved from its stored triangulation, its hyperbolic structure is verified
(verify_hyperbolicity: a certified, positively oriented solution of the gluing equations) and its volume is
computed as a verified interval, which must contain the floating-point volume in the results.
"""
import sys, os, json, csv, glob, time
import snappy
import sage.all  # noqa: F401  (SnapPy's verified computations need Sage)


def verify(tri, vol=None):
    """(ok, lower, upper, bits, error): verified hyperbolicity and volume interval for one triangulation"""
    M = snappy.Manifold(tri)
    err = ''
    for bits in (53, 106, 212):
        try:
            ok = M.verify_hyperbolicity(bits_prec=bits)[0]
            if not ok:
                err = 'not verified'
                continue
            V = M.volume(verified=True, bits_prec=bits)
            lo, hi = float(V.lower()), float(V.upper())
            if vol is not None and not (lo - 1e-9 <= vol <= hi + 1e-9):
                return False, lo, hi, bits, f'stored volume {vol} outside interval'
            return True, lo, hi, bits, ''
        except (KeyboardInterrupt, SystemExit):
            raise
        except BaseException as e:   # cysignals reports a crash inside PARI as SignalError, a BaseException
            err = f'{type(e).__name__}: {str(e)[:100]}'
            M = snappy.Manifold(tri)
            M.randomize()
    return False, None, None, None, err


def task(item):
    t0 = time.time()
    ok, lo, hi, bits, err = verify(item['tri'], item.get('volume'))
    out = {k: v for k, v in item.items() if k != 'tri'}
    out.update(verified=ok, vol_lower=lo, vol_upper=hi, bits=bits, error=err, seconds=round(time.time() - t0, 1))
    return out


def run(items, out, procs):
    # resumable: rows already in `out` are kept and not recomputed.  Workers come from a fork server, since
    # workers forked from the pool's handler thread can crash inside PARI.
    from multiprocessing import get_context
    key = lambda r: (r['word'], str(r['lift_cusp']))
    done = list(csv.DictReader(open(out, newline=''))) if os.path.exists(out) else []
    seen = {key(r) for r in done}
    todo = [it for it in items if key(it) not in seen]
    n_ok = sum(r['verified'] == 'True' for r in done)
    print(f'{len(done)} rows already in {out} ({n_ok} verified), {len(todo)} to do', flush=True)
    fields = list(done[0]) if done else None
    with open(out, 'a', newline='') as fh, get_context('forkserver').Pool(procs, maxtasksperchild=50) as P:
        wr = csv.DictWriter(fh, fieldnames=fields) if fields else None
        for i, r in enumerate(P.imap_unordered(task, todo), start=len(done)):
            if wr is None:
                wr = csv.DictWriter(fh, fieldnames=list(r))
                wr.writeheader()
            wr.writerow(r)
            fh.flush()
            n_ok += bool(r['verified'])
            if not r['verified'] or (i + 1) % 100 == 0:
                print(f"{i + 1}/{len(items)} verified {n_ok}" + ('' if r['verified'] else f"  FAILED {r}"), flush=True)
    print(f'DONE {out}: {n_ok}/{len(items)} verified', flush=True)


def rows_of(workdir):
    for f in sorted(glob.glob(os.path.join(workdir, '*.json'))):
        d = json.load(open(f))
        for r in d.get('rows', []):
            if r.get('volume') not in ('', None) and str(r.get('hyperbolic', True)) == 'True':
                yield d['word'], r


if __name__ == '__main__':
    workdir, out = sys.argv[1], sys.argv[2]
    procs = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    items = [dict(word=w, lift_cusp=r['lift_cusp'], m=r['m'], surface_word=r.get('surface_word'),
                  volume=r['volume'], tri=r['_tri']) for w, r in rows_of(workdir)]
    items.sort(key=lambda x: len(x['tri']))
    print(len(items), 'hyperbolic lifts', flush=True)
    run(items, out, procs)
