#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

# details="+single_atom.echo=9"
# outfile="potential_generator.new_ion"
# for ion in {0..9} ; do
#   echo "# start ion .$ion"
#   $exe +verbosity=7 \
#     -test potential_generator. \
#          +potential_generator.test.ion=.$ion \
#          +potential_generator.max.scf=19 \
#         $details \
#         > $outfile$ion
#   for atom in {0..1} ; do
#   grep 'update_valence_states found a true 3[sp]-eigenstate' $outfile$ion | grep "a#$atom" | tail -2
#     grep ' core     1s   2.0 E=   ' $outfile$ion | grep "a#$atom" | tail -1
#     grep "potential projection for atom #$atom v_00 =" $outfile$ion | head -1
#     grep "potential projection for atom #$atom v_00 =" $outfile$ion | tail -1
#   done
# done

$exe +verbosity=7 \
    -test potential_generator. \
        +electrostatic.solver=load \
        +electrostatic.potential.from.file=v_es.fft.dat \
        +eigensolver=davidson \
        > potential_generator.out.load.dav
