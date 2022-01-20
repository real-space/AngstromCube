#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43


method=C
for Z in {1..1} ; do
  out_Z=data.charge_deficit/$method.z$Z.out

  for energy in {10..10} ; do
    E=`echo $energy \* 0.1 | bc -l`
    out=$method.out
    ./a43 -t single_atom \
            +single_atom.nn.limit=2 \
            +single_atom.test.Z=$Z \
            +single_atom.partial.wave.energy=$E \
            +single_atom.partial.wave.method=$method \
            +single_atom.check.overlap.eigenvalues=-10 \
            +verbosity=7 \
            "$@" > $out

    echo ""                                            >> $out_Z
    echo "# method= $method"                           >> $out_Z
    echo "# energy split= $E Hartree"                  >> $out_Z
    grep 'eigenvalues of the s-overlap operator' $out  >> $out_Z
    grep 'eigenvalues of the p-overlap operator' $out  >> $out_Z
    grep 'eigenvalues of the d-overlap operator' $out  >> $out_Z
    grep 'eigenvalues of the f-overlap operator' $out  >> $out_Z
    echo ""                                            >> $out_Z

  done # energy

done # Z
