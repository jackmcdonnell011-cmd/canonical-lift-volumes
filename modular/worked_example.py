"""The worked example of the paper: the code word w = XXXYXYY followed from start to finish.

usage:  python worked_example.py [word]        (default XXXYXYY)
Opens SnapPy's link viewer for the diagram if it is available (the figure in the paper was drawn there).
"""
import sys

import snappy

from dt_codes import canonical, get_dt, is_primitive
from volumes import calculate_length, word_matrix

w = sys.argv[1] if len(sys.argv) > 1 else 'XXXYXYY'
n = len(w)
print(f'word w = {w}')
print(f'primitive: {is_primitive(w)};  canonical rotation: {canonical(w) == w}')

shifts = [w[i:] + w[:i] for i in range(n)]
rank = {s: r + 1 for r, s in enumerate(sorted(shifts))}
ranks = [rank[s] for s in shifts]
print('ranks of the cyclic shifts S_0..S_{}: {}'.format(n - 1, ', '.join(map(str, ranks))))
print('positions visited on the branch line: ' + ' -> '.join(map(str, ranks + [ranks[0]])))
left = sum(s[0] == 'X' for s in shifts)
print(f'positions 1..{left} on the left branch, {left + 1}..{n} on the right branch')

dt = get_dt(w)
print(f'DT code: {dt}')
L = snappy.Link('DT: ' + dt)
print(f'crossings: {len(L.crossings)}')
M = L.exterior()
print(f'volume: {float(M.volume()):.4f}  ({M.solution_type()})')

(a, b), (c, d) = word_matrix(w)
print(f'M_w = [[{a}, {b}], [{c}, {d}]],  tr = {a + d},  length = 2 arccosh({(a + d) / 2:g}) = {calculate_length(w):.3f}')

try:
    L.view()
except Exception:
    pass
