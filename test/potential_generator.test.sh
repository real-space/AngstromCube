#!/usr/bin/env bash

exe=../src/a43

### simple cubic
printf " 1 \n#cell 12 12 12 i i i\n" > atoms.xyz
echo "Ar   0 0 0" >> atoms.xyz

(cd ../src/ && make -j) && \
$exe +verbosity=11 \
    -test potential_generator. \
        +repeat.eigensolver=1 \
        +eigensolver=cg +conjugate_gradients.max.iter=2 \
        +potential_generator.grid.spacing=0.251 \
        +electrostatic.solver=Bessel0 \
        +element_He="1s* 2 2p | 1.5 sigma .75" \
        +element_Ne="2s 2 2p 6 | 1.8 sigma .7" \
        +element_Ar="3s* 2 3p* 6 3d | 1.6 sigma .8" \
        +element_Xe="5s* 2 5p* 6 5d | 2.24 sigma .62" \
        +occupied.bands=4 \
        +element_B="2s* 2 2p 1 0 3d | 1.2 sigma .7" \
        +element_Al="3s* 2 3p* 1.1 0 3d | 1.8 sigma 1.1 Z= 14" \
        +potential_generator.use.bessel.projection=5 \
        +single_atom.local.potential.method=sinc \
        +single_atom.init.echo=7 \
        +single_atom.init.scf.maxit=1 \
        +single_atom.echo=7 \
        +logder.start=2 +logder.stop=1 \
        +bands.per.atom=10 \
        +potential_generator.max.scf=1 \
        +atomic.valence.decay=0 \
        > potential_generator.out

#         +electrostatic.solver=mg +electrostatic.potential.to.file=v_es.mg.dat \

#         +element_He="1s* 2 2p | 1.5 sigma .75" \
#         +electrostatic.solver=load +electrostatic.potential.from.file=v_es.mg.dat \
#         +electrostatic.solver=load +electrostatic.potential.from.file=v_es.fft.dat \
#         +electrostatic.solver=mg \
#         +element_Ne="2s* 2 2p* 6 3d | 1.8 sigma .7" \


### Results:
### Al-P dimer: Spectrum  -0.636   -0.504   -0.458   -0.219 | -0.216   -0.216   -0.182   -0.134 Ha with #cell 10.5835 10.5835 12.7003 p p p

### He simple cubic at Gamma point: -0.664    0.298    1.092    1.094    1.100    1.173    1.181    2.524    2.528    2.533
### B="2s* 2 2p* 1 0 3d | 1.2 sigma .7" +occupied.bands=1.5    -0.266    0.535    0.538    0.538    0.954    1.276    1.276    1.722    1.756    1.756
### B="2s* 2 2p 1 0 3d | 1.2 sigma .7"  +occupied.bands=1.5    -0.265    0.829    0.846    0.846    0.954    1.285    1.285    1.944    2.009    2.009

### dipole in Si chain experiment:
#
# printf " 4 \n#cell 4.2334 4.2334 8.4668 p p p\n" > atoms.xyz
# echo "Al   0 0 -1.05835" >> atoms.xyz
# echo "P    0 0  1.05835" >> atoms.xyz
# echo "Si   0 0  3.17505" >> atoms.xyz
# echo "Si   0 0 -3.17505" >> atoms.xyz
# 
# (cd ../src/ && make -j) && \
# $exe +verbosity=9 \
#     -test potential_generator. \
#         +eigensolver=none \
#         +potential_generator.grid.spacing=0.251 \
#         +electrostatic.solver=mg +electrostatic.potential.to.file=v_es.mg.dat \
#         +element_He="1s* 2 2p | 1.5 sigma .75" \
#         +element_Ne="2s* 2 2p* 6 3d | 1.8 sigma .7" \
#         +element_Al="3s* 2 3p* 1.1 0 3d | 1.8 sigma 1.1 Z= 14" \
#         +element_Si="3s* 2 3p* 2   0 3d | 1.8 sigma 1.1" \
#          +element_P="3s* 2 3p* 2.9 0 3d | 1.8 sigma 1.1 Z= 14" \
#         +potential_generator.use.bessel.projection=0 \
#         +single_atom.local.potential.method=sinc \
#         +single_atom.init.echo=7 \
#         +single_atom.init.scf.maxit=9 \
#         +single_atom.echo=7 \
#         +logder.start=2 +logder.stop=1 \
#         +bands.per.atom=4 \
#         +potential_generator.max.scf=13 \
#         +atomic.valence.decay=0 \
#         > potential_generator.out
