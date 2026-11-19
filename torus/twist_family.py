"""Volumes along the twist family B^2 A^k on the once-punctured torus.

usage: python -I twist_family.py <kmax> <workdir> <out.csv> [processes]

The curves B^2 A^k (k = 1..kmax, in the free basis a, b of pi_1(S); capitals are inverses) are located on
the modular surface with twist_words.py, which writes their code words and degrees m and tests exactly
whether each one fills S.  The filling ones are run with run_surface.py (all checks), and the volume of each
member is read off the row with the same code word and m.
"""
import csv, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
kmax, workdir, out = int(sys.argv[1]), sys.argv[2], sys.argv[3]
procs = sys.argv[4] if len(sys.argv) > 4 else '4'
py = [sys.executable, '-I']
words, fmap = workdir + '_words.txt', workdir + '_map.csv'
subprocess.run(py + [os.path.join(HERE, 'twist_words.py'), 'torus11', str(kmax), 'BB:A', '--out', words,
                     '--map', fmap], check=True, stdout=subprocess.DEVNULL)
fam = list(csv.DictReader(open(fmap)))
L = max(int(r['code_length']) for r in fam if r['fills'] == 'True')
subprocess.run(py + [os.path.join(HERE, 'run_surface.py'), '--surface', 'torus11', '--L', str(L), '--word-file',
                     words, '--workdir', workdir, '--out', workdir + '_results.csv', '--procs', procs,
                     '--timeout', '10800'], check=True)
rows = []
for r in fam:
    vol = ''
    if r['fills'] == 'True':
        d = json.load(open(os.path.join(workdir, r['modular_word'] + '.json')))
        vol = next((x['volume'] for x in d.get('rows', []) if int(x['m']) == int(r['m'])), '')
    rows.append(dict(k=r['k'], surface_word=r['surface_word'], modular_word=r['modular_word'], m=r['m'],
                     fills=r['fills'], volume=vol))
    print(r['k'], r['surface_word'], r['fills'], vol)
with open(out, 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
    wr.writeheader()
    wr.writerows(rows)
