"""Code words to run, from the census of closed curves on the surface (surface_census.py).

usage: python -I prep_words.py <census.csv> <words.txt> <maxL.txt>

Every filling curve on S is the lift of a modular code word w with degree m; this writes the distinct code
words w of the filling curves and the largest code length m|w| that occurs (the length bound for run_surface).
"""
import csv, sys

census, words_out, maxl_out = sys.argv[1:4]
words, max_l = set(), 0
for r in csv.DictReader(open(census, encoding='utf-8')):
    if r['fills'] == 'True':
        words.add(r['modular_word'])
        max_l = max(max_l, int(r['code_length']))
open(words_out, 'w').write('\n'.join(sorted(words)) + '\n')
open(maxl_out, 'w').write(str(max_l))
print(f'{len(words)} code words, code length <= {max_l}')
