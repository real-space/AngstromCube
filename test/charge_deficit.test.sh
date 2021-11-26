#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

for energy in {1..399} ; do
  E=`echo $energy \* 0.01 | bc -l`
  method=c
  ./a43 -t single_atom \
          +single_atom.nn.limit=2 \
          +single_atom.partial.wave.energy=$E \
          +single_atom.partial.wave.method=$method \
          +verbosity=7 \
          "$@" > out

# for method in c m e r 1 2 ; do
#   ./a43 -t single_atom \
#           +single_atom.nn.limit=1 \
#           +single_atom.partial.wave.method=$method \
#           +single_atom.partial.wave.energy=0.1 \
#           +logder.step=1e-3 \
#           +element_Cu="4s* 1 0 4p* 2e-99 3d 10 | 2.0 sigma .61" \
#           +verbosity=7 \
#           "$@" > out
#         +single_atom.partial.wave.energy=-8 \ ### use energy derivatives instead of excited states

  echo ""
  echo "# method= $method"
  echo "# energy split= $E Hartree"
  grep 'eigenvalues of charge deficit operator for ell=s' out
  grep 'eigenvalues of charge deficit operator for ell=p' out
  grep 'eigenvalues of charge deficit operator for ell=d' out
  grep 'eigenvalues of charge deficit operator for ell=f' out
  echo ""

done
