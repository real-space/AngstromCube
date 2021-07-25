#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

for Z in {12..12}; do

    ## compute the self-consistent solution using atom_core
    reference_file=single_atom.scf.$Z.ref
    $exe +verbosity=3 \
        -t atom_core \
          +atom_core.test.Z=$Z \
         > $reference_file

    ## compute the self-consistent solution using single_atom
    result_file=single_atom.scf.$Z.out
    $exe +verbosity=2 \
        -t single_atom \
          +single_atom.test.Z=$Z \
          +single_atom.nn.limit=1 \
          +single_atom.test.atomic.valence.density=1 \
          +single_atom.relax.partial.waves=0 \
          +single_atom.init.max.scf=66 \
          "$@" > $result_file

done
