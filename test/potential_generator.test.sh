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

$exe +verbosity=17 \
    -test potential_generator. \
        +eigensolver=cg \
        +element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1" \
        +element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma 1.1" \
        +electrostatic.solver=mg \
        -v
#         > potential_generator.out.cg
#         +electrostatic.solver=load +electrostatic.potential.from.file=v_es.mg.dat \

#         +electrostatic.solver=mg +electrostatic.potential.to.file=v_es.mg.dat \
#         +electrostatic.solver=load +electrostatic.potential.from.file=v_es.fft.dat \
