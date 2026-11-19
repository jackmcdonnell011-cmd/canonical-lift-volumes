"""Precompute and cache the Ghys map rho for every modular word with |w| <= N.
Usage: python -I build_rho_cache.py --maxlen 10 --procs 12 [--cache rho_cache]
Resumable: words already in the cache are skipped.  Prints one line per word."""
import os, sys, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import psl2


def words_upto(n):
    out = []
    for L in range(2, n + 1):
        for k in range(2 ** L):
            w = ''.join('XY'[(k >> i) & 1] for i in range(L))
            if 'X' in w and 'Y' in w and w == psl2.canon_word(w) and psl2.primitive_root(w) == w:
                out.append(w)
    return out


def work(args):
    w, cache = args
    t = time.time()
    from lifts import ModularLink
    try:
        ml = ModularLink(w, cache_dir=cache)
        return f'{w} ok {ml.solver} {time.time()-t:.1f}s'
    except Exception as e:
        return f'{w} FAILED {type(e).__name__}: {e} {time.time()-t:.1f}s'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--maxlen', type=int, default=10)
    ap.add_argument('--procs', type=int, default=4)
    ap.add_argument('--hard', action='store_true', help='more attempts + extra solver (for words that FAILED)')
    ap.add_argument('--words', default=None, help='file with one word per line (instead of all words up to --maxlen)')
    ap.add_argument('--cache', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rho_cache'))
    a = ap.parse_args()
    if a.hard:
        os.environ['GN_HARD'] = '1'
    ws = open(a.words).read().split() if a.words else words_upto(a.maxlen)
    todo = [w for w in ws if not os.path.exists(os.path.join(a.cache, w + '.json'))]
    print(f'{len(ws)} words, {len(todo)} to solve', flush=True)
    # longest first so the slow ones start early
    todo.sort(key=lambda w: -len(w))
    from multiprocessing import Pool
    with Pool(a.procs) as P:
        for line in P.imap_unordered(work, [(w, a.cache) for w in todo]):
            print(line, flush=True)
