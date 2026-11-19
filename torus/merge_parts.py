"""Merge the chunk CSVs of one surface and run the cross-row duplicate checks.

  python -I merge_parts.py --surface torus11 --L 22 --word-file words.txt --workdir work --out results.csv

Reads every <workdir>/<word>.json (workdir defaults to work_<surface>_L<L>), so it
does not matter how the chunks were split.  For every pair of hyperbolic rows whose
volumes agree to 1e-8 it runs SnapPy's is_isometric_to on the two single-lift
complements and requires an isometry taking the lifted-knot cusp to the lifted-knot
cusp.  Pairs without one are marked SUSPECT.  Writes the merged CSV and a .summary.txt.
"""
import os
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import argparse
import csv
import json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--surface', required=True)
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--workdir', default=None)
    ap.add_argument('--filling-only', action='store_true', help='the run used run_surface --filling-only')
    ap.add_argument('--word-file', default=None, help='the run used run_surface --word-file')
    a = ap.parse_args()
    import snappy
    import psl2
    from run_surface import word_list
    from pipeline import inverse_word, swap_word, J_map, cusp_iso
    from gn_volumes import SURFACES
    G = SURFACES[a.surface]()
    Jm = J_map(G)
    workdir = a.workdir or os.path.join(HERE, f'work_{a.surface}_L{a.L}')
    words = sorted(set(open(a.word_file).read().split())) if a.word_file else word_list(a.surface, a.L)
    if a.filling_only:
        from run_surface import filling_words
        words = filling_words(a.surface, a.L, words)
    missing, skipped, failed, rows = [], [], [], []
    for w in words:
        fn = os.path.join(workdir, w + '.json')
        if not os.path.exists(fn):
            missing.append(w)
            continue
        res = json.load(open(fn))
        if not res['ok']:
            e = res.get('error', '')
            skipped.append((w, 'rho' if e.startswith('RuntimeError: rho not found') else
                               'timeout' if e.startswith('timeout') else 'other'))
            continue
        if not all(res['checks'].values()):
            failed.append(w)
        rows.extend(res['rows'])
    # cross-row equal-volume checks
    hyp = sorted([r for r in rows if r['hyperbolic']], key=lambda r: r['_vol'])
    cross = []
    cacheF = {}

    def F(r):
        key = (r['word'], r['lift_cusp'])
        if key not in cacheF:
            cacheF[key] = snappy.Manifold(r['_tri'])
        return cacheF[key]
    # cusp-isometry is an equivalence relation, so each row is compared only with one representative of every
    # class found so far in its equal-volume run; the pairs below are then read off the classes.  A test that fails with an error (None) opens a new class, so its
    # pairs are reported as not isometric (SUSPECT), as before.
    cls = {}
    i = 0
    while i < len(hyp):
        j = i
        while j + 1 < len(hyp) and abs(hyp[j + 1]['_vol'] - hyp[i]['_vol']) <= 1e-8:
            j += 1
        reps = []
        for k in range(i, j + 1):
            R = hyp[k]
            for c, Q in enumerate(reps):
                if cusp_iso(F(Q), Q['_Kpos'], F(R), R['_Kpos']) is True:
                    cls[k] = (i, c)
                    break
            else:
                cls[k] = (i, len(reps))
                reps.append(R)
        i = j + 1
    for i in range(len(hyp)):
        for j in range(i + 1, len(hyp)):
            A, B = hyp[i], hyp[j]
            if abs(A['_vol'] - B['_vol']) > 1e-8:
                break
            iso = cls[i] == cls[j]
            wa, wb = A['word'], B['word']
            if wa == wb:
                rel = 'same word, different component'
            elif inverse_word(wa) == wb:
                rel = 'inversion (kappa vs kappa^-1)'
            elif swap_word(wa) == wb:
                rel = 'J (mirror)' + ('' if Jm is not None else ', Gamma not J-invariant')
            elif psl2.canon_word(wa[::-1]) == wb:
                rel = 'J o inversion' + ('' if Jm is not None else ', Gamma not J-invariant')
            else:
                rel = 'no known symmetry of S (e.g. a mapping class)'
            cross.append((wa, A['lift_cusp'], wb, B['lift_cusp'], round(A['_vol'], 9), iso, rel))
            for R, other in ((A, f'{wb}#{B["lift_cusp"]}'), (B, f'{wa}#{A["lift_cusp"]}')):
                t = f'equal vol to {other}: cusp-isometry={iso}, relation={rel}'
                if iso is not True:
                    t = 'SUSPECT ' + t
                R['status'] = t if R['status'] == 'ok' else R['status'] + '; ' + t
    cols = ['word', 'modular_word_length', 'm', 'surface_word_length', 'hyperbolic_length_on_S', 'volume',
            'multiplicity', 'gamma_trace', 'checks', 'gamma_element', 'hyperbolic', 'solution_type', 'status', 'lift_cusp', 'surface_word']
    rows.sort(key=lambda r: (r['surface_word_length'], r['hyperbolic_length_on_S'], r['word'], r['lift_cusp']))
    with open(a.out, 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        wr.writeheader()
        wr.writerows(rows)
    clean = [r for r in rows if r['hyperbolic'] and 'SUSPECT' not in r['status'] and 'EXCLUDED' not in r['status']]
    suspect = [(r['word'], r['lift_cusp']) for r in rows if 'SUSPECT' in r['status']]
    summ = (f'surface={a.surface} L={a.L} Gamma={G.name} index={G.d}\n'
            f'words={len(words)} missing(not yet run)={len(missing)} {missing}\n'
            f'skipped={len(skipped)} ' + ' '.join(
                f'{k}={[w for w, r in skipped if r == k]}' for k in ('rho', 'timeout', 'other')) + '\n'
            f'failed_checks={len(failed)} {failed}\n'
            f'rows={len(rows)} hyperbolic_rows={sum(1 for r in rows if r["hyperbolic"])} '
            f'clean_hyperbolic_rows={len(clean)} distinct_volumes={len({round(r["_vol"], 7) for r in clean})}\n'
            f'suspect={len(suspect)} {suspect}\n'
            f'equal_volume_pairs={len(cross)}\n' + ''.join(f'  {c}\n' for c in cross))
    open(os.path.splitext(a.out)[0] + '.summary.txt', 'w').write(summ)
    print(summ)


if __name__ == '__main__':
    main()
