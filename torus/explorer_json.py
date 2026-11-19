"""Explorer entry for a surface with length measured in the surface's own generators.

usage: python -I explorer_json.py results_N11.csv results_N11.summary.txt census_N11.csv 1,1 "Once-punctured torus (via the commutator subgroup Γ′)" ../explorer/torus11.json

census.csv (surface_census.py) lists every closed curve up to some surface word length, with whether
it fills.  results.csv (run_torus.sh) has the volumes, with the surface word of each lift.  Each
filling curve is matched to its row by surface word; for a normal Gamma one row stands for all lifts
of a word (equal volumes, checked by the pipeline), so it is matched by (modular word, m) instead.
Rows whose complements are isometric (merge_parts.py) get the same manifold-type number.
"""
import ast, csv, json, sys
csv.field_size_limit(2**31 - 1)      # status fields can list many equal-volume partners

res_path, summary_path, census_path, key, label, out = sys.argv[1:7]

rows = [r for r in csv.DictReader(open(res_path, encoding='utf-8'))
        if r['volume'] and str(r['hyperbolic']) == 'True'
        and (all(c.endswith('pass') for c in r['checks'].split())
             or r['checks'] == 'via the reversed curve')      # reverse_rho.py: checks passed on the reversed word
        and 'EXCLUDED' not in r['status']]
by_sw = {r['surface_word']: r for r in rows if r.get('surface_word')}
by_wm = {(r['word'], r['m']): r for r in rows if int(r['multiplicity']) > 1}

parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
for line in open(summary_path, encoding='utf-8'):
    line = line.strip()
    if line.startswith("('"):
        w1, i1, w2, i2, _vol, iso, _rel = ast.literal_eval(line)
        if iso:
            parent[find((w1, i1))] = find((w2, i2))

census = [c for c in csv.DictReader(open(census_path, encoding='utf-8')) if c.get('fills') == 'True']
pts, missing = [], []
for c in census:
    r = by_sw.get(c['surface_word']) or by_wm.get((c['modular_word'], c['m']))
    if r is None:
        missing.append(c['surface_word'])
        continue
    pts.append((c, r))
# a hyperbolic row must be a filling curve
fill_words = {c['surface_word'] for c in census}
N = max(int(c['surface_word_length']) for c in census)
unexpected = [r['surface_word'] for r in rows if r.get('surface_word') and len(r['surface_word']) <= N
              and r['surface_word'] not in fill_words]

roots = {}
for c, r in sorted(pts, key=lambda p: (float(p[1]['volume']), p[1]['word'])):
    roots.setdefault(find((r['word'], int(r['lift_cusp'] or 0))), len(roots) + 1)
s = {
    'label': label,
    'maxL': N,
    'fullUpTo': N,
    'counts': {},
    'w': [c['surface_word'] for c, r in pts],
    'code': [r['word'] + (f'^{r["m"]}' if int(r['m']) > 1 else '') for c, r in pts],
    'wl': [int(c['surface_word_length']) for c, r in pts],
    'len': [round(float(r['hyperbolic_length_on_S']), 3) for c, r in pts],
    'vol': [round(float(r['volume']), 4) for c, r in pts],
    'cls': [roots[find((r['word'], int(r['lift_cusp'] or 0)))] for c, r in pts],
    'families': {},
    'wlName': 'word length in the generators of π₁(S)',
    'note': ('Each point is one closed curve on the surface, written as a cyclic word in a free basis of '
             'π₁(S) (the side pairings of a fundamental polygon made of Farey triangles; letters A, B, C, … '
             'and their inverses a, b, c, …), counted once up to orientation. Only curves whose canonical '
             'lift complement is hyperbolic are shown; these are exactly the curves that fill the surface. '
             'The hover also gives the modular code word of the curve.'),
}
for k in s['wl']:
    s['counts'][k] = s['counts'].get(k, 0) + 1
json.dump({key: s}, open(out, 'w', encoding='utf-8'), ensure_ascii=False)
print(f'{key}: {len(pts)} curves, {len(roots)} manifold types, counts {dict(sorted(s["counts"].items()))}; '
      f'filling curves without a volume: {len(missing)} {missing[:8]}; hyperbolic rows not in the census '
      f'filling list: {len(unexpected)} {unexpected[:8]}')
