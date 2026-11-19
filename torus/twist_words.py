"""Twist families u t^k on a surface: their modular code words, so they can be run with run_surface.py.

usage: python -I twist_families.py <surface> <kmax> <u:t> [<u:t> ...] [--find N] [--out words.txt] [--map map.csv] [--spiral]

Words are written as in the census (lowercase = basis generator, uppercase = inverse).  t must not be
peripheral, so that u t^k is an immersed twist along the closed geodesic t rather than a spiral into a cusp;
--spiral allows a peripheral t, for families of the first kind.
--find N lists candidate families: u t and u t^2 both fill, |u t^2| <= N, t of length 1 or 2 and not peripheral.
"""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
args = sys.argv[1:]
spiral = '--spiral' in args
if spiral:
    args.remove('--spiral')
surface, kmax = args[0], int(args[1])
opt = {}
for key in ('--find', '--out', '--map'):
    if key in args:
        i = args.index(key)
        opt[key] = args[i + 1]
        del args[i:i + 2]
fams = [a.split(':') for a in args[2:]]
sys.argv = [sys.argv[0], surface, '0']
import surface_census as sc
from surface_words import cyclic_reduce, letters, is_primitive_cyclic


def parse(w):
    return [ord(ch) - 96 if ch.islower() else -(ord(ch) - 64) for ch in w]


def info(word):
    u = cyclic_reduce(parse(word))
    if not is_primitive_cyclic(u):
        return None
    loc = sc.locate(u)
    if loc is None:
        return None
    w, m, c, o = loc
    return w, m, c, sc.fe.fills(w, c, m)


if '--find' in opt:
    N = int(opt['--find'])
    r = len(sc.B.gens)
    alpha = [letters([i]) for i in range(1, r + 1)] + [letters([-i]) for i in range(1, r + 1)]
    ts = [a for a in alpha] + [a + b for a in alpha for b in alpha if a != b.swapcase() and a != b]
    ts = [t for t in ts if info(t) is not None]
    us = [a + b for a in alpha for b in alpha if a != b.swapcase()]
    us += [x + c for x in us for c in alpha if x[-1] != c.swapcase()]
    for t in ts:
        for u in us:
            if len(u) + 2 * len(t) > N or u[-1] == t[0].swapcase() or u[0] == t[-1].swapcase() or u[-1] == t[-1] and len(t) == 1:
                continue
            a, b = info(u + t), info(u + t * 2)
            if a and b and a[3] and b[3]:
                print(f'{u}:{t}  code lengths {a[1]*len(a[0])}, {b[1]*len(b[0])}')
    sys.exit()

words, rows = set(), []
for u, t in fams:
    assert spiral or info(t) is not None, f'{t} is peripheral (use --spiral)'
    for k in range(1, kmax + 1):
        x = info(u + t * k)
        if x is None:
            continue
        w, m, c, f = x
        rows.append(dict(family=f'{u}:{t}', k=k, surface_word=u + t * k, modular_word=w, m=m, code_length=m * len(w), fills=f))
        if f:
            words.add(w)
        print(f'{u}{t}^{k}: code length {m * len(w)}, m={m}, fills={f}')
if '--out' in opt:
    open(opt['--out'], 'w').write('\n'.join(sorted(words, key=len)) + '\n')
if '--map' in opt:
    with open(opt['--map'], 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
