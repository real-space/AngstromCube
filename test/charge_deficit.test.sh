#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

for energy in {1..399} ; do

  E=`echo $energy \* 0.01 | bc -l`
  ./a43 -t single_atom +single_atom.nn.limit=2 -V +single_atom.partial.wave.energy=$E > out
  grep 'charge deficit operator for ell=s' out
  grep 'charge deficit operator for ell=p' out
  grep 'charge deficit operator for ell=d' out
  grep 'charge deficit operator for ell=f' out

done
