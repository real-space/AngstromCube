#!/usr/bin/env bash

exe=../src/a43


### Al atom
# geometry_file=atoms.xyz
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "Al   0 0 0" >> $geometry_file
# project_base=pg_new.Al-atom

### C-sc
# geometry_file=atoms.xyz
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "C  0 0 0" >> $geometry_file
# project_base=pg_new.C-sc

### C-dimer:
geometry_file=atoms.xyz
printf " 2 \n#cell 8 8 8 p p p \n" > $geometry_file
echo "C  0 0 -0.65" >> $geometry_file
echo "C  0 0  0.65" >> $geometry_file
project_base=pg_new.C-dimer
### QE-spectrum          -19.6629 -10.5856  -7.2179  -7.2179  -7.1693  -0.4972  -0.4972  -0.1544 eV (self-consistent result)
### exe using davidson   -22.522  -11.837   -3.817   -0.825   -0.341    1.797    1.913    2.057 2.414 eV (not self-consistent)
### exe using CG         -22.694  -12.207   -4.013   -4.013   -2.165    0.242    1.671    1.707   1.721   2.005 2.005 2.195 2.195 3.065 3.363 4.106 4.137 4.165 4.180 4.329
### exe using SHO4       -21.2913 -9.80567 -7.60901 -7.60899 -7.39258 -0.223547 -0.223528 11.133 eV
### exe using 10Ry       -22.6197 -12.1007 -4.01218 -4.01207 -2.19573 0.234374 1.69084 1.7138    133.364 133.388 eV

### H-sc
# geometry_file=atoms.xyz
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "H  0 0 0" >> $geometry_file
# project_base=pg_new.H-sc


### Al-P dimer
# geometry_file=atoms.xyz
# printf " 2 \n#cell 4.233418 4.233418 8.466836 p p p \n" > $geometry_file
# printf " 2 \n#cell 10.5835 10.5835 12.7003 p p p \n" > $geometry_file
# printf " 2 \n#cell 21.16708996 21.16708996 25.400507952 p p p \n" > $geometry_file
# echo "Al   0 0 -1.058354498" >> $geometry_file
# echo "P    0 0  1.058354498" >> $geometry_file
# project=potential_generator.AlP.pw.out

### Al fcc bulk
# geometry_file=atoms.xyz
# printf " 4 \n#cell 4.0 4.0 4.0 p p p \n" > $geometry_file
# echo "Al   -1.0 -1.0 -1.0" >> $geometry_file
# echo "Al    1.0  1.0 -1.0" >> $geometry_file
# echo "Al    1.0 -1.0  1.0" >> $geometry_file
# echo "Al   -1.0  1.0  1.0" >> $geometry_file
# project_base=pg.Al-fcc

### P sc bulk
# geometry_file=atoms.xyz
# printf " 1 \n#cell 3.0 3.0 3.0 p p p \n" > $geometry_file
# echo "P 0 0 0" >> $geometry_file
# project_base=potential_generator.P-sc

### Al sc bulk
# geometry_file=atoms.xyz
# printf " 1 \n#cell 3.0 3.0 3.0 p p p \n" > $geometry_file
# echo "Al 0 0 0" >> $geometry_file
# project_base=potential_generator.Al-sc

# geometry_file=C_chain.xyz
# printf " 1 \n#cell 1.420282 10 10 p p p \n" > $geometry_file
# echo "C   0 0 0" >> $geometry_file
# project_base=potential_generator.$geometry_file.sho.out

# geometry_file=graphene.xyz
# project_base=potential_generator.graphene.sho.out

nkpoints=2
scale_k=1
scale_p=1
scale_h=1
scale_s=1

### generate a control file
cat > control.sh << EOF
# verbosity of log output
verbosity=7

# where to find the coordinates
geometry.file=$geometry_file

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


# configuration for basis=grid
bands.per.atom=10
repeat.eigensolver=15
eigensolver=davidson
# this option is not active
grid_hamiltonian.floating.point.bits=32

# for the CG method
conjugate_gradients.max.iter=19

# for start wave functions use SHO functions with larger sigma spread
start.waves.scale.sigma=9
atomic.valence.decay=0

# show energies in units of electronVolt
output.energy.unit=eV

# configuration for basis=sho
sho_hamiltonian.test.sigma=1.0
sho_hamiltonian.test.sigma.asymmetry=1
sho_hamiltonian.test.kpoints=$nkpoints
sho_hamiltonian.scale.kinetic=$scale_k
sho_hamiltonian.scale.potential=$scale_p
sho_hamiltonian.scale.nonlocal.h=$scale_h
sho_hamiltonian.scale.nonlocal.s=$scale_s
sho_hamiltonian.floating.point.bits=32


dense_solver.test.overlap.eigvals=1

# configuration for basis=pw
pw_hamiltonian.test.kpoints=$nkpoints
pw_hamiltonian.scale.kinetic=$scale_k
pw_hamiltonian.scale.potential=$scale_p
pw_hamiltonian.scale.nonlocal.h=$scale_h
pw_hamiltonian.scale.nonlocal.s=$scale_s
pw_hamiltonian.floating.point.bits=32

#! maybe we could unify some sho_hamiltonian and pw_hamiltonian options
EOF


for spacing in `seq 4 1 4`; do
  project=$project_base.grid$spacing
  echo "# start calculation $project"

  (cd ../src/ && make -j) && \
  $exe -test potential_generator. \
        +control.file=control.sh \
        +basis=grid \
        > $project.out
        ./spectrum.sh $project.out > $project.spectrum.dat
done

for numax in `seq 4 1 4`; do
  project=$project_base.sho$numax
  echo "# start calculation $project"

  (cd ../src/ && make -j) && \
  $exe -test potential_generator. \
        +control.file=control.sh \
        +basis=sho \
        +sho_hamiltonian.test.numax=$numax \
        > $project.out
        ./spectrum.sh $project.out > $project.spectrum.dat
done

for ecut in `seq 1 1 5`; do
  project=$project_base.pw$ecut
  echo "# start calculation $project"

  (cd ../src/ && make -j) && \
  $exe -test potential_generator. \
        +control.file=control.sh \
        +basis=pw \
        +pw_hamiltonian.cutoff.energy=$ecut \
        > $project.out
        ./spectrum.sh $project.out > $project.spectrum.dat
done
