#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43


# for i in {0..9}; do
#   echo $i
#   grep 'lowest eigenvalues for ell=' single_atom.out.energy$i
# done
# exit

# $exe -t single_atom. -V \
#          +single_atom.partial.wave.method=e \
#          +single_atom.partial.wave.energy.derivative=0 \
#          +single_atom.partial.wave.energy=1.0

### stable (Cu)
# $exe -t single_atom. -V \
#          +single_atom.partial.wave.method=m \
#          +single_atom.partial.wave.energy.derivative=0 \
#          +single_atom.partial.wave.energy=1.0 \
#          > single_atom.out


# for Z in {13..15}; do
#   $exe +verbosity=7 \
#     -test single_atom. \
#         +single_atom.test.Z=$Z \
#         +single_atom.partial.wave.method=m \
#         +single_atom.partial.wave.energy.derivative=0 \
#         +single_atom.partial.wave.energy=1.0 \
#         > single_atom.out.$Z.new
# done

# for Z in {58..70}; do
#   $exe +verbosity=7 \
#     -test single_atom. \
#         +single_atom.test.Z=$Z \
#         +single_atom.partial.wave.method=m \
#         +single_atom.partial.wave.energy.derivative=0 \
#         +single_atom.partial.wave.energy=1.0 \
#         > single_atom.out.$Z.new
# done


# out=single_atom.out.warn
# echo -n "# " > $out
# date     >> $out
# echo " " >> $out
# echo " " >> $out
# echo " " >> $out
# 
# for Z in {1..118}; do
#   $exe +verbosity=7 \
#     -test single_atom. \
#         +single_atom.test.Z=$Z \
#         +single_atom.partial.wave.method=m \
#         +single_atom.partial.wave.energy.derivative=0 \
#         +single_atom.partial.wave.energy=1.0 \
#         +single_atom.nn.limit=1 \
#         +logder.start=2 +logder.stop=-3 \
#         >> $out
# done

# Z=29
# out="single_atom.out.$Z.try"
# echo -n "# " > $out
# date     >> $out
# echo " " >> $out
# 
#   $exe +verbosity=7 \
#     -test single_atom. \
#         +single_atom.test.Z=$Z \
#         +single_atom.partial.wave.method=m \
#         +single_atom.test.atomic.valence.density=1 \
#         +single_atom.partial.wave.energy.derivative=0 \
#         +single_atom.nn.limit=1 \
#         +logder.start=3 +logder.stop=-2 \
#         +single_atom.init.scf.maxit=9 \
#         >> $out
# exit

# Z=29
# out="single_atom.out.$Z.logder"
# echo -n "# " > $out
# date     >> $out
# echo " " >> $out
# echo " " >> $out
# echo " " >> $out
# 
# for E in $(seq -2 0.01 1); do
#   echo "Energy parameter = $E"
#   $exe +verbosity=7 \
#     -test single_atom. \
#         +single_atom.test.Z=$Z \
#         +single_atom.partial.wave.method=m \
#         +single_atom.test.spherical.valence=1 \
#         +single_atom.partial.wave.energy.derivative=0 \
#         +single_atom.partial.wave.energy.parameter=$E \
#         +logder.start=$E +logder.stop=$E \
#         +single_atom.nn.limit=1 \
#         >> $out
# done
# eval
# grep '#  emm-averaged  0 hamiltonian  ' $out | awk '{print $5 }'  > single_atom.out.emm0
# grep '#  emm-averaged  2 hamiltonian  ' $out | awk '{print $7 }'  > single_atom.out.emm1
# grep '#  emm-averaged  4 hamiltonian  ' $out | awk '{print $9 }'  > single_atom.out.emm2
# grep '#  emm-averaged  5 hamiltonian  ' $out | awk '{print $10 }' > single_atom.out.emm3
# 
# grep '#  spherical 0s hamiltonian  ' $out | awk '{print $5 }'   > single_atom.out.sph0
# grep '#  spherical 0p hamiltonian  ' $out | awk '{print $7 }'   > single_atom.out.sph1
# grep '#  spherical 0d hamiltonian  ' $out | awk '{print $9 }'   > single_atom.out.sph2
# grep '#  spherical 0f hamiltonian  ' $out | awk '{print $10 }'  > single_atom.out.sph3
# 
# grep '#  logarithmic_derivative check ' $out -A1 > single_atom.out.logder 


for Z in {1..10}; do
  ## a separate out file for each Z
  out_Z=single_atom.out.$Z.try
  echo -n "# using sigma_config    " > $out_Z
  date >> $out_Z

  echo "# Z = $Z" ## show progress in terminal
  $exe +verbosity=17 \
    -test single_atom. \
        +smooth.radial.grid.from=0 \
        +single_atom.test.Z=$Z \
        +single_atom.partial.wave.energy.derivative=0 \
        +single_atom.partial.wave.energy=0 \
        +single_atom.partial.wave.method=e \
        +single_atom.test.atomic.valence.density=1 \
        +single_atom.from.sigma.config=1 \
        +single_atom.nn.limit=2 \
        +logder.start=-2 +logder.stop=-3 +logder.step=.001 \
        >> $out_Z
done
### seems like the eigenvalues of the charge deficit for ell=0 are always -1 and some positive value up to Z <= 86 and for Z=119,120



exit




### We find that when we use a fixed energy parameter, 
###     the eigenenergies found in eigenstate_analysis match almost exactly 
###     that energy parameter as long as be are bound (E < 0)

out=single_atom.out.conf
echo -n "# using sigma_config    " > $out
date     >> $out
echo " " >> $out

for Z in {117..117}; do
  ## either make a separate out file for each Z
  out_Z=single_atom.out.$Z.conf.new
  echo -n "# using sigma_config    " > $out_Z
  date >> $out_Z
  ## or concat all into the global out file
  # out_Z=$out

  echo "# Z = $Z" ## show progress in terminal
  $exe +verbosity=7 \
    -test single_atom. \
        +single_atom.test.Z=$Z \
        +single_atom.partial.wave.method=m \
        +single_atom.partial.wave.energy.derivative=0 \
        +single_atom.partial.wave.energy=1.0 \
        +single_atom.test.atomic.valence.density=1 \
        +single_atom.from.sigma.config=1 \
        +logder.start=-2 +logder.stop=1 +logder.step=.001 \
        +single_atom.nn.limit=1 \
        +show.state.diagram=-1 \
        >> $out_Z
done
## We find stability execept for ell=2 of Z=7--9 and ell=0 of Z=80,84--92,117--120
##   which we need to analyze further for the "well-configured" elements (Z<=86)

