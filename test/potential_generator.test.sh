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

$exe +verbosity=9 \
    -test potential_generator. \
        +eigensolver=cg +conjugate_gradients.max.iter=2 \
        +repeat.eigensolver=9 \
        +electrostatic.solver=load +electrostatic.potential.from.file=v_es.mg.dat \
        +element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma 1.1" \
         +element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1" \
        +potential_generator.use.bessel.projection=0 \
        +single_atom.local.potential.method=sinc \
        +single_atom.init.echo=6 \
        +single_atom.echo=6 \
        +logder.start=2 +logder.stop=1 \
        +bands.per.atom=4 \
      > potential_generator.out

#         +element_He="1s* 2 2p | 1.5 sigma .75" \
#         +electrostatic.solver=mg +electrostatic.potential.to.file=v_es.mg.dat \
#         +electrostatic.solver=load +electrostatic.potential.from.file=v_es.fft.dat \
#         +electrostatic.solver=mg \


### Results:
### Al-P dimer: Spectrum  -0.636   -0.504   -0.458   -0.219 | -0.216   -0.216   -0.182   -0.134 Ha
