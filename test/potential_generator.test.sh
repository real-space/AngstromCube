#!/usr/bin/env bash

exe=../src/a43
geometry_file=atoms.xyz


# project_base=pg.Al-atom
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "Al   0 0 0" >> $geometry_file

# project_base=pg.C-sc
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "C  0 0 0" >> $geometry_file

project_base=pg.C-atom
printf " 1 \n#cell 6 6 6 p p p \n" > $geometry_file
echo "C  0 0 0" >> $geometry_file

# project_base=pg.Mg-atom
# printf " 1 \n#cell 6 6 6 p p p \n" > $geometry_file
# echo "Mg  0 0 0" >> $geometry_file


# project_base=pg.C-dimer
# printf " 2 \n#cell 8 8 8 p p p \n" > $geometry_file
# echo "C  0 0 -0.65" >> $geometry_file
# echo "C  0 0  0.65" >> $geometry_file
### QE-spectrum          -19.6629 -10.5856  -7.2179  -7.2179  -7.1693  -0.4972  -0.4972  -0.1544 eV (self-consistent result)
### exe using davidson   -22.522  -11.837   -3.817   -0.825   -0.341    1.797    1.913    2.057 2.414 eV (not self-consistent)
### exe using CG         -22.694  -12.207   -4.013   -4.013   -2.165    0.242    1.671    1.707   1.721   2.005 2.005 2.195 2.195 3.065 3.363 4.106 4.137 4.165 4.180 4.329
### exe using SHO4       -21.2913 -9.80567 -7.60901 -7.60899 -7.39258 -0.223547 -0.223528 11.133 eV
### exe using 10Ry       -22.6197 -12.1007 -4.01218 -4.01207 -2.19573 0.234374 1.69084 1.7138    133.364 133.388 eV

# project_base=pg.H-sc
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "H  0 0 0" >> $geometry_file


### Al-P dimer
# printf " 2 \n#cell 4.233418 4.233418 8.466836 p p p \n" > $geometry_file
# printf " 2 \n#cell 10.5835 10.5835 12.7003 p p p \n" > $geometry_file
# printf " 2 \n#cell 21.16708996 21.16708996 25.400507952 p p p \n" > $geometry_file
# echo "Al   0 0 -1.058354498" >> $geometry_file
# echo "P    0 0  1.058354498" >> $geometry_file
# project_base=potential_generator.AlP

### Al fcc bulk
# printf " 4 \n#cell 4.0 4.0 4.0 p p p \n" > $geometry_file
# echo "Al   -1.0 -1.0 -1.0" >> $geometry_file
# echo "Al    1.0  1.0 -1.0" >> $geometry_file
# echo "Al    1.0 -1.0  1.0" >> $geometry_file
# echo "Al   -1.0  1.0  1.0" >> $geometry_file
# project_base=pg.Al-fcc

### P sc bulk
# printf " 1 \n#cell 3.0 3.0 3.0 p p p \n" > $geometry_file
# echo "P 0 0 0" >> $geometry_file
# project_base=potential_generator.P-sc

### Al sc bulk
# printf " 1 \n#cell 3.0 3.0 3.0 p p p \n" > $geometry_file
# echo "Al 0 0 0" >> $geometry_file
# project_base=potential_generator.Al-sc

### C_chain
# printf " 1 \n#cell 1.420282 10 10 p p p \n" > $geometry_file
# echo "C   0 0 0" >> $geometry_file
# project_base=potential_generator.C_chain.sho

### generate a control file
cat > control.sh << EOF

# verbosity of the log files
verbosity=7

# show energies in units of electronVolt
output.energy.unit=eV

# grid spacing of the dense grid
potential_generator.grid.spacing=0.23622

# max number of self-consistency iterations
potential_generator.max.scf=1

# Poisson solver
electrostatic.solver=fft


# configuration of atomic PAW setups
#element_H="1s 1 0 | 0.9 sigma .41"
#element_C="2s 2 2p 2 0 | 1.2 sigma .38 numax 2 V=sinc"
#element_C="2s 2 2p 2 0 | 1.2 numax 1 sigma .445 V=parabola"
element_C="2s 2 2p 2 0 | 1.2 numax 2 sigma .38 V=sinc"
#element_C="2s* 2 2p 2 0 | 1.2 numax 2 sigma .38 V=sinc"
#element_C="2s** 2 2p* 2 0 3d 4f | 1.2 numax 4 sigma .314327 V=sinc"
#element_C="2s*** 2 2p** 2 0 3d** 4f* 5g* 6h 7i | 1.2 numax 6 sigma .33055552 V=sinc"
#element_C="2s 2 2p 2 0 | 1.2 numax 6 sigma .33055552 V=sinc"
#element_C="2s 2 2p 2 0 | 1.2 sigma .445 numax 9 V=parabola"
#element_Mg="2s 2 3s 2 2p 6 3p 2e-99 3d | 1.96 numax 3 sigma .60636 V=sinc"
#element_Mg="2s 2 3s* 2 2p 6 3p 2e-99 3d | 1.96 numax 4 sigma .578 V=sinc"
#element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma .5 V=parabola"
#element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1 V=sinc"
#single_atom.local.potential.method=sinc
single_atom.nn.limit=4
single_atom.partial.wave.method=energy_ordering
single_atom.echo=7
single_atom.init.echo=7
single_atom.optimize.sigma=0
single_atom.init.scf.maxit=0
# logarithmic derivatives
logder.unit=Ha
logder.start=-2
logder.stop=1
#logder.step=1e-3

# configuration for basis=grid
bands.per.atom=10
# DEVEL option
devel.occupied.bands=8
eigensolver=cg
repeat.eigensolver=3
conjugate_gradients.max.iter=19
# for start wave functions use SHO functions with larger sigma spread
start.waves.scale.sigma=9
atomic.valence.decay=0
export.waves=waves.dat
#start.waves=waves.dat

# configuration for basis=sho
# spread of the SHO basis functions in Bohr
sho_hamiltonian.test.sigma=1.0

# configuration for basis=sho or basis=pw
hamiltonian.test.kpoints=2
hamiltonian.scale.kinetic=1
hamiltonian.scale.potential=1
hamiltonian.scale.nonlocal.h=1
hamiltonian.scale.nonlocal.s=1
hamiltonian.floating.point.bits=32

# configuration for basis=pw
pw_hamiltonian.solver=direct
pw_hamiltonian.density=1
#pw_hamiltonian.solver=iterative
#pw_hamiltonian.solver=both
davidson_solver.max.iterations=9

# also compute the eigenvalues of the overlap matrix?
dense_solver.test.overlap.eigvals=0

EOF


for spacing in `seq 1 1 0`; do
  project=$project_base.grid$spacing
  (cd ../src/ && make -j) && \
  echo "# start calculation $project" && \
  $exe -test potential_generator. \
        +control.file=control.sh \
        +basis=grid \
        > $project.out
        ./spectrum.sh $project.out > $project.spectrum.dat
#         +element_C="2s 2 2p 2 0 | 1.2 sigma .445 numax $spacing V=parabola" \
done

for numax in `seq 4 2 3`; do
  project=$project_base.sho$numax
  (cd ../src/ && make -j) && \
  echo "# start calculation $project" && \
  $exe -test potential_generator. \
        +control.file=control.sh \
        +basis=sho \
        +sho_hamiltonian.test.numax=$numax \
        > $project.out
        ./spectrum.sh $project.out > $project.spectrum.dat
done

for ecut in `seq 2 1 2`; do
  project=$project_base.pw$ecut
  (cd ../src/ && make -j) && \
  echo "# start calculation $project" && \
  $exe -test potential_generator. \
        +control.file=control.sh \
        +basis=pw \
        +pw_hamiltonian.cutoff.energy=$ecut \
        > $project.out
        ./spectrum.sh $project.out > $project.spectrum.dat
done
