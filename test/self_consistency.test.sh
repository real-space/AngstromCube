#!/usr/bin/env bash

exe=../src/a43
geometry_file=atoms.xyz

# project_base=scf.C-sc
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "C 0 0 0" >> $geometry_file

# project_base=scf.vacuum
# printf " 1 \n#cell 4 4 4 i i i \n" > $geometry_file
# echo "__  0 0 0" >> $geometry_file

project_base=scf.C-atom
printf " 1 \n#cell 8 8 8 i i i \n" > $geometry_file
echo "C  0 0 0" >> $geometry_file

# project_base=scf.C-dimer
# printf " 2 \n#cell 8 8 8 i i i \n" > $geometry_file
# echo "C  -0.65 0 0" >> $geometry_file
# echo "C   0.65 0 0" >> $geometry_file
## test translational invariance
# echo "C  0 0 -0.525" >> $geometry_file
# echo "C  0 0  0.775" >> $geometry_file

# project_base=potential_generator.AlP
# printf " 2 \n#cell 4.233418 4.233418 8.466836 p p p \n" > $geometry_file
# printf " 2 \n#cell 10.5835 10.5835 12.7003 p p p \n" > $geometry_file
# printf " 2 \n#cell 21.16708996 21.16708996 25.400507952 p p p \n" > $geometry_file
# echo "Al   0 0 -1.058354498" >> $geometry_file
# echo "P    0 0  1.058354498" >> $geometry_file

# project_base=scf.Al-fcc
# printf " 4 \n#cell 4.0 4.0 4.0 p p p \n" > $geometry_file
# echo "Al   -1.0 -1.0 -1.0" >> $geometry_file
# echo "Al    1.0  1.0 -1.0" >> $geometry_file
# echo "Al    1.0 -1.0  1.0" >> $geometry_file
# echo "Al   -1.0  1.0  1.0" >> $geometry_file

### Cu LDA lattice constant from PHYSICAL REVIEW B 79, 085104 􏰀(2009􏰁), al. et Blaha
# project_base=scf.Cu-fcc
# printf " 4 \n#cell 3.522 3.522 3.522 p p p \n" > $geometry_file
# echo "Cu   -.8805 -.8805 -.8805" >> $geometry_file
# echo "Cu    .8805  .8805 -.8805" >> $geometry_file
# echo "Cu    .8805 -.8805  .8805" >> $geometry_file
# echo "Cu   -.8805  .8805  .8805" >> $geometry_file

### diamond LDA lattice constant 3.536 Ang from PHYSICAL REVIEW B 79, 085104 􏰀(2009􏰁), al. et Blaha
# project_base=scf.C-diamond
# printf " 8 \n#cell 3.536 3.536 3.536 p p p \n" > $geometry_file
# ## 3.536 / 8 == 0.442, 0.442 * 3 == 1.326
# echo "C  -1.326 -1.326 -1.326" >> $geometry_file
# echo "C   0.442  0.442 -1.326" >> $geometry_file
# echo "C  -0.442 -0.442 -0.442" >> $geometry_file
# echo "C   1.326  1.326 -0.442" >> $geometry_file
# echo "C  -1.326  0.442  0.442" >> $geometry_file
# echo "C   0.442 -1.326  0.442" >> $geometry_file
# echo "C  -0.442  1.326  1.326" >> $geometry_file
# echo "C   1.326 -0.442  1.326" >> $geometry_file

# project_base=potential_generator.C_chain.sho
# printf " 1 \n#cell 1.420282 10 10 p p p \n" > $geometry_file
# echo "C   0 0 0" >> $geometry_file



### generate a control file
cat > control.sh << EOF
################################################################
#
#   a43 control file
#
################################################################
# use # for silent comments and #! for comments in the log file

## general verbosity of the log file
verbosity=7

## display energies in custom units {Ha, Ry, eV}
output.energy.unit=eV
## display distances in custom units {Bohr, nm, Ang}
output.length.unit=Bohr

## grid spacing or number of grid points of the dense grid
grid.spacing.unit=Ang
#grid.spacing=0.083334
#grid.spacing=0.10416667
grid.points=24

## max. number of self-consistency iterations
self_consistency.max.scf=3
self_consistency.mix.density=0.25
atomic.valence.decay=1
## compute the Fermi level {exact, linearized}
fermi.level=exact

## analyze the potentials up to vtot (DEBUG)
#potential_generator.use.bessel.projection=0
#potential_generator.use.direct.projection=0

## Poisson equation solver {multigrid, none, cg, sd, fft}
electrostatic.solver=fft


##
## configuration of atomic PAW setups
##

## configurations of elements
#element_C="2s* 2 2p 2 0 | 1.2 numax 2 sigma .4304 V=parabola"
element_C="2s* 2 2p* 2 0 | 1.2 numax 3 sigma .345353 V=parabola"

#element___="1s 1e-6 2p 0 | 1.0 numax 1 sigma .5 V=parabola"

# configuration of atomic PAW setups
#element_H="1s 1 0 | 0.9 sigma .308 V=parabola"
#element_He="1s 2 2p | 1.5 numax 1 sigma .537535"
#element_Li="2s 1 0 2p 2e-99 | 2.0 numax 1 sigma .7752 V=parabola"
#element_Li="2s 1 0 2p 2e-99 | 2.0 numax 2 sigma .612475 V=sinc"
#element_Li="2s 1 0 2p 2e-99 | 2.0 numax 1 sigma .8088 V=sinc"
#element_C="2s 2 2p 2 0 | 1.2 numax 1 sigma .445 V=parabola"
#element_C="2s 2 2p 2 0 | 1.2 sigma .38 numax 2 V=sinc"
#element_C="2s* 2 2p 2 0 | 1.2 numax 2 sigma .38 V=sinc"
#single_atom.partial.wave.energy=1.1 good for carbon with 2s*
#element_C="2s** 2 2p* 2 0 3d 4f | 1.2 numax 4 sigma .314327 V=sinc"
#element_C="2s*** 2 2p** 2 0 3d** 4f* 5g* 6h 7i | 1.2 numax 6 sigma .33055552 V=sinc"
#element_C="2s 2 2p 2 3d 4f 5g 6h 7i | 1.2 numax 15 sigma .197208 V=sinc"
#element_C="2s 2 2p 2 | 1.2 numax 15 sigma .19416 V=parabola"
#element_C="2s 2 2p 2 0 3d 4f 5g 6h 7i | 1.2 numax 6 sigma .33055552 V=sinc"
#element_C="2s 2 2p 2 0 | 1.2 numax 6 sigma .33055552 V=sinc"
#element_C="2s 2 2p 2 0 | 1.2 sigma .445 numax 9 V=parabola"
#element_Mg="2s 2 3s 2 2p 6 3p 2e-99 3d | 1.96 numax 3 sigma .60636 V=sinc"
#element_Mg="2s 2 3s* 2 2p 6 3p 2e-99 3d | 1.96 numax 4 sigma .578 V=sinc"
element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma .5 V=parabola"
element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1 V=sinc"

### Copper: for the occupation setting 4s 1 3d 10, the d-states should be 1.5 eV below the s-state
### element_Cu="4s 2 4p 2e-99 3d 5 4 | 2.0 sigma .651 V=sinc" (seems to work well with GRID, but
###                          fails with PW, s-state below d-states, ss-density matrix entry tiny)
### element_Cu="4s 1 0 4p 2e-99 3d 10 | 2.2 numax 2 sigma .71857 V=sinc" --> s-state below d-states, as well

## partial wave method {energy_ordering, recreate_second, classical, ...}
single_atom.partial.wave.method=recreate_second

## relax partial wave in every SCF iteration {1:yes, 0:no}
single_atom.relax.partial.waves=0

## special verbosity for PAW setup
single_atom.init.echo=5

## special verbosity for PAW update
single_atom.echo=3

## default for local potential method {parabola, sinc}
#single_atom.local.potential.method=sinc

## bit mask for the first 50 atoms, -1:all, 1:only atom#0, 5:atoms#0 and #2 but not #1, ...
single_atom.echo.mask=1

## optimize sigma for the occupied projectors {0:no, 1:yes, -1:optimize and show but don't use it, -10:optimize and exit}
single_atom.optimize.sigma=1

## debug options
#single_atom.synthetic.density.matrix=1
single_atom.init.max.scf=0

## export PAW data in paw_xml format for GPAW or ABINIT
single_atom.export.xml=0

###! the same radial grid is used to represent true and smooth quantities
### smooth.radial.grid.from=0

## to check logarithmic derivatives ensure start < stop
logder.unit=eV
logder.start=90
logder.stop=50
logder.step=1e-2

##
## general configurations for the Kohn-Sham equation
##

## number of Kohn-Sham states per atom in the system
bands.per.atom=10

## sampling of the Brillouin zone
#hamiltonian.kmesh.echo=9
#hamiltonian.kmesh.x=3
#hamiltonian.kmesh=0
## hamiltonian.floating.point.bits=64
hamiltonian.floating.point.bits=32

## configuration for basis=grid
# method of the grid eigensolver {cg, Davidson, none, explicit}
#grid.eigensolver=none
#grid.eigensolver=explicit
grid.eigensolver=cg
#conjugate_gradients.max.iter=1
#grid.eigensolver.repeat=1

## for start wave functions use SHO functions with larger sigma spread
start.waves.scale.sigma=6

## load start waves from file, store wave functions to file
#start.waves=$project_base.waves.dat
store.waves=$project_base.waves.dat

## configuration for basis=pw
#plane_wave.solver {auto, both, direct, iterative}
plane_wave.solver=direct

#plane_wave.solver=iterative
#plane_wave.iterative.solver.ratio=8.0
#plane_wave.iterative.solver=cg
#conjugate_gradients.max.iter=2
#plane_wave.max.cg.iterations=12
#davidson_solver.max.iterations=7

## also compute the eigenvalues of the overlap matrix?
#dense_solver.test.overlap.eigvals=0

# Export the Hamiltonian {0:not, 1:yes, -1:yes+abort} in format {xml, json}
export.hamiltonian=1
export.hamiltonian.format=xml

# Make structure_solver produce the same as potential_generator
# structure_solver.complex=1

control.show=1

EOF


for spacing in `seq 2 1 2`; do
  project=$project_base.grid$spacing
  (cd ../src/ && make -j) && \
  echo "# start calculation $project" && \
  $exe -test self_consistency \
        +control.file=control.sh \
        +basis=grid \
        "$@" > $project.out
  ./spectrum.sh $project.out > $project.spectrum.dat
done

for numax in `seq 4 2 0`; do
  project=$project_base.sho$numax
  (cd ../src/ && make -j) && \
  echo "# start calculation $project" && \
  $exe -test self_consistency \
        +control.file=control.sh \
        +basis=sho \
        +sho_hamiltonian.test.numax=$numax \
        +sho_hamiltonian.test.sigma=.5 \
        "$@" > $project.out
  ./spectrum.sh $project.out > $project.spectrum.dat
done

for ecut in `seq 5 2 0`; do
  project=$project_base.pw$ecut
  (cd ../src/ && make -j) && \
  echo "# start calculation $project" && \
  $exe -test self_consistency \
        +control.file=control.sh \
        +basis=pw \
        +plane_wave.cutoff.energy=$ecut \
        "$@" > $project.out
  ./spectrum.sh $project.out > $project.spectrum.dat
done

