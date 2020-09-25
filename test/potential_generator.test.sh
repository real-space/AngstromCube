#!/usr/bin/env bash

exe=../src/a43

nkpoints=2

### Al atom
# geometry_file=atoms.xyz
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "Al   0 0 0" >> $geometry_file
# out_file_base=pg_new.Al-atom

### C-sc
# geometry_file=atoms.xyz
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "C  0 0 0" >> $geometry_file
# out_file_base=pg_new.C-sc

### C-dimer:
geometry_file=atoms.xyz
printf " 2 \n#cell 8 8 8 p p p \n" > $geometry_file
echo "C  0 0 -0.65" >> $geometry_file
echo "C  0 0  0.65" >> $geometry_file
out_file_base=pg_new.C-dimer
### QE-spectrum          -19.6629 -10.5856  -7.2179  -7.2179  -7.1693  -0.4972  -0.4972  -0.1544 eV (self-consistent result)
### exe using davidson   -22.522  -11.837   -3.817   -0.825   -0.341    1.797    1.913    2.057 2.414 eV (not self-consistent)
### exe using CG         -22.694  -12.207   -4.013   -4.013   -2.165    0.242    1.671    1.707   1.721   2.005 2.005 2.195 2.195 3.065 3.363 4.106 4.137 4.165 4.180 4.329
### exe using SHO4       -21.2913 -9.80567 -7.60901 -7.60899 -7.39258 -0.223547 -0.223528 11.133 eV
### exe using 10Ry       -22.6197 -12.1007 -4.01218 -4.01207 -2.19573 0.234374 1.69084 1.7138    133.364 133.388 eV

### H-sc
# geometry_file=atoms.xyz
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "H  0 0 0" >> $geometry_file
# out_file_base=pg_new.H-sc


### Al-P dimer
# geometry_file=atoms.xyz
# printf " 2 \n#cell 4.233418 4.233418 8.466836 p p p \n" > $geometry_file
# printf " 2 \n#cell 10.5835 10.5835 12.7003 p p p \n" > $geometry_file
# printf " 2 \n#cell 21.16708996 21.16708996 25.400507952 p p p \n" > $geometry_file
# echo "Al   0 0 -1.058354498" >> $geometry_file
# echo "P    0 0  1.058354498" >> $geometry_file
# out_file=potential_generator.AlP.pw.out

### Al fcc bulk
# geometry_file=atoms.xyz
# printf " 4 \n#cell 4.0 4.0 4.0 p p p \n" > $geometry_file
# echo "Al   -1.0 -1.0 -1.0" >> $geometry_file
# echo "Al    1.0  1.0 -1.0" >> $geometry_file
# echo "Al    1.0 -1.0  1.0" >> $geometry_file
# echo "Al   -1.0  1.0  1.0" >> $geometry_file
# out_file_base=pg.Al-fcc

### P sc bulk
# geometry_file=atoms.xyz
# printf " 1 \n#cell 3.0 3.0 3.0 p p p \n" > $geometry_file
# echo "P 0 0 0" >> $geometry_file
# out_file_base=potential_generator.P-sc

### Al sc bulk
# geometry_file=atoms.xyz
# printf " 1 \n#cell 3.0 3.0 3.0 p p p \n" > $geometry_file
# echo "Al 0 0 0" >> $geometry_file
# out_file_base=potential_generator.Al-sc

# geometry_file=C_chain.xyz
# printf " 1 \n#cell 1.420282 10 10 p p p \n" > $geometry_file
# echo "C   0 0 0" >> $geometry_file
# out_file_base=potential_generator.$geometry_file.sho.out

# geometry_file=graphene.xyz
# out_file_base=potential_generator.graphene.sho.out

### generate a control file
cat > control.sh << EOF
# grid spacing of the dense grid
potential_generator.grid.spacing=0.251

# max number of self-consistency iterations
potential_generator.max.scf=1

# Poisson solver
electrostatic.solver=fft

# DEVEL option
occupied.bands=4

# configuration of atomic PAW setups
element_H="1s 1 0 | 0.9 sigma .41"
element_C="2s 2 2p 2 0 | 1.2 sigma .5"
element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma .5"
element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1"
single_atom.local.potential.method=sinc
single_atom.init.echo=7
single_atom.init.scf.maxit=1
single_atom.echo=1
# logarithmic derivatives
logder.unit=Ha
logder.start=2
logder.stop=1

# this option is not active
grid_hamiltonian.floating.point.bits=32

# iterative eigenvalue solver options
bands.per.atom=10
repeat.eigensolver=15
eigensolver=davidson

# for the CG method
conjugate_gradients.max.iter=19

# for start wave functions use SHO functions with larger sigma spread
start.waves.scale.sigma=9
atomic.valence.decay=0

# show energies in units of electronVolt
output.energy.unit=eV
EOF

scale_k=1
scale_p=1
scale_h=1
scale_s=1

for spacing in `seq 4 1 4`; do
  out_file=$out_file_base.grid$spacing.out
  echo "# start calculation $out_file"

  (cd ../src/ && make -j) && \
  $exe +verbosity=7 \
    -test potential_generator. \
        +geometry.file=$geometry_file \
        +control.file=control.sh \
        +basis=grid \
        > $out_file
        ./spectrum.sh $out_file > $out_file.spectrum.dat
done
exit
#         +potential_generator.grid.spacing=0.251 \
#         +electrostatic.solver=fft \
#         +occupied.bands=4 \
#         +element_H="1s 1 0 | 0.9 sigma .41" \
#         +element_C="2s 2 2p 2 0 | 1.2 sigma .5" \
#         +element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma .5" \
#          +element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1" \
#         +single_atom.local.potential.method=sinc \
#         +single_atom.init.echo=7 \
#         +single_atom.init.scf.maxit=1 \
#         +single_atom.echo=1 \
#         +logder.start=2 +logder.stop=1 \
#         +bands.per.atom=10 \
#         +potential_generator.max.scf=1 \
#         +grid_hamiltonian.floating.point.bits=32 \
#         +repeat.eigensolver=15 \
#         +eigensolver=davidson \
#         +conjugate_gradients.max.iter=19 \
#         +start.waves.scale.sigma=9 \
#         +atomic.valence.decay=0 \
#         +output.energy.unit=eV \

for numax in `seq 4 1 4`; do
  out_file=$out_file_base.sho$numax.out
  echo "# start calculation $out_file"

  (cd ../src/ && make -j) && \
  $exe +verbosity=7 \
    -test potential_generator. \
        +geometry.file=$geometry_file \
        +potential_generator.grid.spacing=0.251 \
        +electrostatic.solver=fft \
        +occupied.bands=4 \
        +element_H="1s 1 0 | 0.9 sigma .41" \
        +element_C="2s 2 2p 2 0 | 1.2 sigma .5" \
        +element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma .5" \
         +element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1" \
        +single_atom.local.potential.method=sinc \
        +single_atom.init.echo=7 \
        +single_atom.init.scf.maxit=1 \
        +single_atom.echo=1 \
        +logder.start=2 +logder.stop=1 \
        +bands.per.atom=4 \
        +potential_generator.max.scf=1 \
        +basis=sho \
        +sho_hamiltonian.test.numax=$numax \
        +sho_hamiltonian.test.sigma=1.0 \
        +sho_hamiltonian.test.sigma.asymmetry=1 \
        +sho_hamiltonian.test.kpoints=$nkpoints \
        +sho_hamiltonian.scale.kinetic=$scale_k \
        +sho_hamiltonian.scale.potential=$scale_p \
        +sho_hamiltonian.scale.nonlocal.h=$scale_h \
        +sho_hamiltonian.scale.nonlocal.s=$scale_s \
        +sho_hamiltonian.floating.point.bits=32 \
        +dense_solver.test.overlap.eigvals=1 \
        +output.energy.unit=eV \
        > $out_file
        ./spectrum.sh $out_file > $out_file.spectrum.dat
done

for ecut in `seq 1 1 5`; do
  out_file=$out_file_base.pw$ecut.out
  echo "# start calculation $out_file"

  (cd ../src/ && make -j) && \
  $exe +verbosity=7 \
    -test potential_generator. \
        +geometry.file=$geometry_file \
        +potential_generator.grid.spacing=0.251 \
        +electrostatic.solver=fft \
        +occupied.bands=4 \
        +element_H="1s 1 0 | 0.9 sigma .41" \
        +element_C="2s 2 2p 2 0 | 1.2 sigma .5" \
        +element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma .5" \
         +element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1" \
        +single_atom.local.potential.method=sinc \
        +single_atom.init.echo=7 \
        +single_atom.init.scf.maxit=1 \
        +single_atom.echo=1 \
        +logder.start=2 +logder.stop=1 \
        +bands.per.atom=4 \
        +potential_generator.max.scf=1 \
        +basis=pw \
        +pw_hamiltonian.cutoff.energy=$ecut \
        +pw_hamiltonian.test.kpoints=$nkpoints \
        +pw_hamiltonian.scale.kinetic=$scale_k \
        +pw_hamiltonian.scale.potential=$scale_p \
        +pw_hamiltonian.scale.nonlocal.h=$scale_h \
        +pw_hamiltonian.scale.nonlocal.s=$scale_s \
        +pw_hamiltonian.floating.point.bits=32 \
        +dense_solver.test.overlap.eigvals=1 \
        +output.energy.unit=eV \
        > $out_file
        ./spectrum.sh $out_file > $out_file.spectrum.dat
done
exit


#         +electrostatic.solver=load +electrostatic.potential.from.file=v_es.mg.dat \
#         +element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma 1.1" \
#          +element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1" \
#         +electrostatic.solver=mg \
#         +electrostatic.potential.to.file=v_es.mg.dat \
#         +repeat.eigensolver=0 \
#         +eigensolver=cg +conjugate_gradients.max.iter=1 \
#         +start.waves.scale.sigma=5 \
#         +atomic.valence.decay=0 \

#         +sho_hamiltonian.test.green.function=200 \
#         +sho_hamiltonian.test.green.lehmann=200 \

### simple cubic
printf " 1 \n#cell 12 12 12 i i i\n" > atoms.xyz
echo "Xe   0 0 0" >> atoms.xyz

(cd ../src/ && make -j) && \
$exe +verbosity=11 \
    -test potential_generator. \
        +repeat.eigensolver=5 \
        +eigensolver=cg +conjugate_gradients.max.iter=12 \
        +potential_generator.grid.spacing=0.251 \
        +electrostatic.solver=Bessel0 \
        +element_He="1s* 2 2p | 1.5 sigma .75" \
        +element_Ne="2s 2 2p 6 | 1.8 sigma .7" \
        +element_Ar="3s* 2 3p* 6 3d | 1.6 sigma .8" \
        +element_Ar="3s 2 3p 6 | 1.6 sigma .8" \
        +element_Xe="5s* 2 5p* 6 5d | 2.24 sigma .62" \
        +element_Xe="5s 2 5p 6 | 2.24 sigma .62" \
        +occupied.bands=4 \
        +element_B="2s* 2 2p 1 0 3d | 1.2 sigma .7" \
        +element_Al="3s* 2 3p* 1.1 0 3d | 1.8 sigma 1.1 Z= 14" \
        +potential_generator.use.bessel.projection=5 \
        +single_atom.local.potential.method=parabola \
        +single_atom.init.echo=7 \
        +single_atom.init.scf.maxit=1 \
        +single_atom.echo=7 \
        +logder.start=2 +logder.stop=1 \
        +bands.per.atom=4 \
        +potential_generator.max.scf=5 \
        +start.waves.scale.sigma=9 \
        +atomic.valence.decay=0 \
        > potential_generator.out

#         +electrostatic.solver=mg +electrostatic.potential.to.file=v_es.mg.dat \

#         +element_He="1s* 2 2p | 1.5 sigma .75" \
#         +electrostatic.solver=load +electrostatic.potential.from.file=v_es.fft.dat \
#         +electrostatic.solver=mg \
#         +element_Ne="2s* 2 2p* 6 3d | 1.8 sigma .7" \


### Results:
### Al-P dimer: Spectrum  -0.636   -0.504   -0.458   -0.219 | -0.216   -0.216   -0.182   -0.134 Ha with #cell 10.5835 10.5835 12.7003 p p p

### He simple cubic at Gamma point: -0.664    0.298    1.092    1.094    1.100    1.173    1.181    2.524    2.528    2.533
### B="2s* 2 2p* 1 0 3d | 1.2 sigma .7" +occupied.bands=1.5    -0.266    0.535    0.538    0.538    0.954    1.276    1.276    1.722    1.756    1.756
### B="2s* 2 2p 1 0 3d | 1.2 sigma .7"  +occupied.bands=1.5    -0.265    0.829    0.846    0.846    0.954    1.285    1.285    1.944    2.009    2.009
### Ar="3s 2 3p 6 | 1.6 sigma .8"                              -0.642   -0.140   -0.140   -0.140 (seems too high by .25 Ha)
### Xe="5s 2 5p 6 | 2.24 sigma .62"                            -0.555   -0.234   -0.234   -0.234

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
