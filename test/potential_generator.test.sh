#!/usr/bin/env bash

exe=../src/a43
(cd ../src/ && make -j) && \
$exe +verbosity=9 \
    -test potential_generator. \
        +eigensolver=none \
        +potential_generator.grid.spacing=0.251 \
        +electrostatic.solver=mg +electrostatic.potential.to.file=v_es.mg.dat \
        +element_He="1s* 2 2p | 1.5 sigma .75" \
        +element_Ne="2s* 2 2p* 6 3d | 1.8 sigma .7" \
        +element_Al="3s* 2 3p* 1.5 0 3d | 1.8 sigma 1.1 Z= 14" \
         +element_P="3s* 2 3p* 2.5 0 3d | 1.8 sigma 1.1 Z= 14" \
        +potential_generator.use.bessel.projection=0 \
        +single_atom.local.potential.method=sinc \
        +single_atom.init.echo=7 \
        +single_atom.init.scf.maxit=9 \
        +single_atom.echo=7 \
        +logder.start=2 +logder.stop=1 \
        +bands.per.atom=4 \
        +potential_generator.max.scf=13 \
        +atomic.valence.decay=0 \
        > potential_generator.out

#         +repeat.eigensolver=9 \
#         +eigensolver=cg +conjugate_gradients.max.iter=12 \
#         +element_He="1s* 2 2p | 1.5 sigma .75" \
#         +electrostatic.solver=load +electrostatic.potential.from.file=v_es.mg.dat \
#         +electrostatic.solver=load +electrostatic.potential.from.file=v_es.fft.dat \
#         +electrostatic.solver=mg \


### Results:
### Al-P dimer: Spectrum  -0.636   -0.504   -0.458   -0.219 | -0.216   -0.216   -0.182   -0.134 Ha with #cell 10.5835 10.5835 12.7003 p p p
