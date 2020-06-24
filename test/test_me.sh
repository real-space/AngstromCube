#!/usr/bin/env bash

### this is a run script that makes sure that the binary is up to date with the sources before running the executable
exe=../src/a43
echo "#" $exe "$@"

Z=2

echo -n "# " ## comment out "make: Nothing to be done for `all'."
(cd ../src/ && make -j) && \
  $exe +verbosity=7  \
  -test single_atom. \
       +single_atom.test.Z=$Z \
       +single_atom.test.atomic.valence.density=1 \
       +single_atom.nn.limit=2 \
       +logder.start=-2.5 +logder.step=.001 +logder.stop=1.5 \
       +single_atom.init.scf.maxit=1 \
       +single_atom.optimize.sigma=-1 \
       +element_H="1s* 1 0 2p 2e-99 | 0.9 sigma .623" \
       +element_He="1s* 2 2p | 1.5 sigma .75" \
       +element_B="2s* 2 2p* 1 0 3d | 1.2 sigma .7" \
       +element_C="2s* 2 2p* 2 0 3d | 1.2 sigma .8" \
       +element_N="2s* 2 2p* 3 0 3d | 1.0 sigma .8" \
       +element_O="2s* 2 2p* 3 1 3d | 1.13 sigma .9" \
       +element_F="2s* 2 2p* 3 2 3d | 1.2 sigma .8" \
       +element_Ne="2s* 2 2p* 6 3d | 1.8 sigma .7" \
       +element_Mg="2s 2 3s 2 2p 6 3p 2e-99 3d | 1.96 sigma .41" \
       +element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1" \
       +element_S="3s* 2 3p* 3 1 3d | 1.7 sigma 1." \
       +element_Ar="3s* 2 3p* 6 3d | 1.6 sigma .8" \
       +element_K="3s 2 4s 1 0 3p* 6 3d | 1.77 sigma .8" \
       +element_Sc="3s 2 4s 2 3p 6 4p 2e-99 3d* 1 0 | 2.32 sigma .6" \
       +element_Cu="4s* 1 1 4p* 2e-99 3d* 9 4f 2e-99 | 2.0 sigma .61" \
       +element_Ge="4s* 2 4p* 2 0 4d | 1.9 sigma 1." \
       +element_As="4s* 2 4p* 3 0 4d | 2.0 sigma 1.1" \
       +element_Se="4s* 2 4p* 3 1 4d | 1.6 sigma 1." \
       +element_Kr="4s* 2 4p* 6 4d | 2.2 sigma .9" \
       +element_Rb="4s 2 5s 1 0 4p* 6 4d | 2.3 sigma 1." \
       +element_u0="7s 2 8s 2 7p 6 8p 2e-99 6d 10 6f 2e-99 | 3. sigma .9" \
       +single_atom.local.potential.method=sinc \
      > single_atom.$Z.sinc

##       +logder.stop=-9 \
## +logder.start=-1 +logder.step=.0001 +logder.stop=0
#        +logder.start=-2.5 +logder.step=.0001 +logder.stop=.5 \
#        +single_atom.local.potential.method=sinc \
