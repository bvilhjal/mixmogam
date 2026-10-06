#!/bin/zsh
# Same-host check of exact LOCO's worker RSS: e4089d8 snapshot against the live 2.0.0.dev6 tree.
W=${0:a:h}
PY=/Users/au507860/anaconda3/envs/ldpred3-accelerate/bin/python
OLD=/Users/au507860/REPOS/mixmogam/benchmarks/results/20261005-phensim-kvik-e4089d8/source/kvik_simulation.py
NEW=/Users/au507860/REPOS/mixmogam/benchmarks/ldak_kvik_comparison.py
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1
one() {  # label driver method
  local cache=$W/cache-$1
  mkdir -p $cache
  ( cd $W/panel/unstructured_null && NUMBA_CACHE_DIR=$cache /usr/bin/time -l $PY $2 --worker $W/panel/unstructured_null --method $3 --seed 20262403 ) 2> $W/time.txt > /dev/null
  echo "$1 $3 $(awk '/real/ {printf "%s", $1} /maximum resident/ {printf " %.1f", $1/1048576}' $W/time.txt) load $(sysctl -n vm.loadavg)"
}
# The e4089d8 driver calls HRATT by its earlier name.
for pair in "exact exact" "kvik hratt"; do
  old_m=${pair%% *}; new_m=${pair##* }
  one e4089d8 $OLD $old_m > /dev/null; one dev6 $NEW $new_m > /dev/null   # warm the Numba caches
  for label in e4089d8 dev6 dev6 e4089d8 e4089d8 dev6 dev6 e4089d8; do
    [[ $label == e4089d8 ]] && one $label $OLD $old_m || one $label $NEW $new_m
  done
done
