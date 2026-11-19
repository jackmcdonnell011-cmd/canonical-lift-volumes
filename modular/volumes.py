"""Hyperbolic volume and length of every closed geodesic on the modular surface up to a given word length.

For each code word w the DT code is passed to SnapPy, whose exterior is the canonical lift complement
M_w = UT(Mod) minus the lift of w = S^3 minus (trefoil and Lorenz knot of w), and volume() solves the gluing
equations.  The length is 2 arccosh(|tr M_w|/2), with the trace computed exactly over the integers.

usage:  python volumes.py <max_length> <out.csv> [processes]
output: word,word_length,trace,hyperbolic_length,volume,solution_type
"""
import csv
import math
import sys
from multiprocessing import Pool

import snappy

from dt_codes import code_words, get_dt

X = ((1, 1), (0, 1))
Y = ((1, 0), (1, 1))


def word_matrix(w):
    """M_w, the product of the matrices of the letters of w"""
    a, b, c, d = 1, 0, 0, 1
    for ch in w:
        (p, q), (r, s) = X if ch == 'X' else Y
        a, b, c, d = a * p + b * r, a * q + b * s, c * p + d * r, c * q + d * s
    return (a, b), (c, d)


def trace(w):
    (a, _), (_, d) = word_matrix(w)
    return abs(a + d)


def calculate_length(word):
    """2 arccosh(|tr M_w| / 2); 0 for parabolic words"""
    t = trace(word)
    return 0.0 if t <= 2 else 2 * math.acosh(t / 2)


def calculate_knot_volume(dt_code_string):
    """volume of the exterior of the link with this DT code, and SnapPy's solution type (the triangulation is
    randomised until the solution has all tetrahedra positively oriented, at most 20 times)"""
    M = snappy.Link('DT: ' + dt_code_string).exterior()
    for _ in range(20):
        if M.solution_type() == 'all tetrahedra positively oriented':
            break
        M.randomize()
    return float(M.volume()), M.solution_type()


def row(w):
    vol, sol = calculate_knot_volume(get_dt(w))
    return [w, len(w), trace(w), repr(calculate_length(w)), repr(vol), sol]


if __name__ == '__main__':
    max_length, out = int(sys.argv[1]), sys.argv[2]
    procs = int(sys.argv[3]) if len(sys.argv) > 3 else None
    words = list(code_words(max_length))
    with Pool(procs) as pool, open(out, 'w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['word', 'word_length', 'trace', 'hyperbolic_length', 'volume', 'solution_type'])
        for r in pool.imap(row, words, chunksize=200):
            wr.writerow(r)
    print(f'{len(words)} words written to {out}')
