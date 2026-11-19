#!/bin/bash
# Canonical lift volumes of every filling closed curve on the once-punctured torus up to surface word length N.
# usage (from this folder):  bash run_torus.sh <N> [processes]      e.g.  bash run_torus.sh 11 8
# output: census_N<N>.csv, results_N<N>.csv (+ .summary.txt), log in run_N<N>.log
#
# 1. surface_census.py lists every closed curve of word length <= N and tests exactly whether it fills.
# 2. run_surface.py builds M_gamma for every lift of every code word with all five checks (resumable: one
#    JSON per word in work_N<N>/).
# 3. Words whose Ghys map rho was not found are retried with more attempts; any still without rho, and lifts
#    that fill but have no hyperbolic structure from SnapPy, take their manifold from the reversed curve
#    (reverse_rho.py, harden.py, reverse_fix.py).
# 4. merge_parts.py merges the words and tests every pair of equal volume for a cusp-preserving isometry;
#    filling_exact.py compares SnapPy's hyperbolicity with the exact filling test on every row.
cd "$(dirname "$0")"
set -e
N=$1; P=${2:-4}; s=torus11
W=work_N$N; LOG=run_N$N.log
python -I surface_census.py $s $N census_N$N.csv | tee -a $LOG
python -I prep_words.py census_N$N.csv words_N$N.txt maxL_N$N.txt | tee -a $LOG
L=$(cat maxL_N$N.txt)
python -I run_surface.py --surface $s --L $L --word-file words_N$N.txt --workdir $W --out results_N$N.csv \
    --procs $P --timeout 3600 >> $LOG 2>&1
python -I merge_parts.py --surface $s --L $L --word-file words_N$N.txt --workdir $W --out results_N$N.csv >> $LOG 2>&1
rho=$(grep -o "rho=\[[^]]*\]" results_N$N.summary.txt | tr -d "[]',=" | sed 's/^rho//')
if [ -n "$rho" ]; then
    printf '%s\n' $rho > rho_retry_N$N.txt
    for w in $rho; do
        python -I -c "import sys; sys.path.insert(0, '.'); from lifts import ModularLink; ModularLink('$w', tries=400, cache_dir='rho_cache')" >> $LOG 2>&1 || true
        [ -f $W/$w.json ] && mv $W/$w.json $W/$w.json.norho
    done
    python -I run_surface.py --surface $s --L $L --word-file rho_retry_N$N.txt --workdir $W \
        --out results_N${N}_rho.csv --procs 2 --timeout 3600 >> $LOG 2>&1
fi
python -I reverse_rho.py $s $W $L 2 >> $LOG 2>&1
python -I harden.py $s $W $P >> $LOG 2>&1
python -I reverse_fix.py $s $W $L $P >> $LOG 2>&1
python -I merge_parts.py --surface $s --L $L --word-file words_N$N.txt --workdir $W --out results_N$N.csv >> $LOG 2>&1
python -I filling_exact.py $s results_N$N.csv | tee -a $LOG
sed -n 2,6p results_N$N.summary.txt
