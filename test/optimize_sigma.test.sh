#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43


method=c
for Z in {1..1} ; do
  out_Z=z$Z.out
  rm -f $out_Z

  for isigma in {20..50} ; do
    sigma=`echo $isigma \* 0.01 | bc -l`
    out=$method.out
    ./a43 -t single_atom \
            +single_atom.test.Z=$Z \
            +single_atom.nn.limit=2 \
            +single_atom.partial.wave.energy=1.0 \
            +single_atom.partial.wave.method=$method \
            +single_atom.check.overlap.eigenvalues=-10 \
            +single_atom.optimize.sigma=0 \
            +element_H="1s* 1 0 2p | 0.9 sigma $sigma numax 15" \
            +verbosity=7 \
            > $out
          

    echo -n "# sigma= $sigma Bohr"                              >> $out_Z
    grep 'orthogonalized coefficients of the 1s-projector' $out >> $out_Z

  done # sigma

done # Z
