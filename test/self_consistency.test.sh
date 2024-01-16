#!/usr/bin/env bash

exe=./a43

# base=scf.H-fcc
# a=2.12132
# c=0.53033
# echo " 4" > $base.xyz
# echo "#cell $a $a $a p p p" >> $base.xyz
# echo "H   -$c -$c -$c" >> $base.xyz
# echo "H    $c  $c -$c" >> $base.xyz
# echo "H    $c -$c  $c" >> $base.xyz
# echo "H   -$c  $c  $c" >> $base.xyz

# base=scf.C-sc
# echo " 1" > $base.xyz
# echo "#cell 3.0 3.0 3.0 p p p" >> $base.xyz
# echo "C 0 0 0" >> $base.xyz

# base=scf.vacuum
# echo " 1" > $base.xyz
# echo "#cell 4 4 4 i i i" >> $base.xyz
# echo "__  0 0 0" >> $base.xyz

base=scf.C-atom
echo " 1" > $base.xyz
echo "#cell 32 32 32 i i i" >> $base.xyz
echo "C  0 0 0" >> $base.xyz

# base=scf.methane
# echo " 5" > $base.xyz
# echo "#cell 16 16 16 i i i" >> $base.xyz
# base=scf.mini-methane
# echo " 5" > $base.xyz
# echo "#cell 8 8 8 i i i" >> $base.xyz
# dist=0.63
# echo "C  0 0 0"                 >> $base.xyz
# echo "H   -$dist -$dist -$dist" >> $base.xyz
# echo "H    $dist  $dist -$dist" >> $base.xyz
# echo "H    $dist -$dist  $dist" >> $base.xyz
# echo "H   -$dist  $dist  $dist" >> $base.xyz

# base=scf.C-chain
# echo " 1" > $base.xyz
# echo "#cell 2 8 8 p i i\n" >> $base.xyz
# echo "C  0 0 0" >> $base.xyz

# base=scf.C-dimer
# echo " 2" > $base.xyz
# echo "#cell 8 8 8 i i i" >> $base.xyz
# echo "C  -0.65 0 0" >> $base.xyz
# echo "C   0.65 0 0" >> $base.xyz
## test translational invariance
# echo "C  0 0 -0.525" >> $base.xyz
# echo "C  0 0  0.775" >> $base.xyz

# base=scf.AlP
# echo " 2" > $base.xyz
# echo "#cell 4.233418 4.233418 8.466836 p p p" >> $base.xyz
# echo "#cell 10.5835 10.5835 12.7003 p p p" >> $base.xyz
# echo "#cell 21.16708996 21.16708996 25.400507952 p p p" >> $base.xyz
# echo "#cell 21.167 21.167 21.167 i i i" >> $base.xyz
# echo "#cell 8.0 8.0 8.0 i i i" >> $base.xyz
# echo "Al   0 0 -1.058354499" >> $base.xyz
# echo "P    0 0  1.058354499" >> $base.xyz

# base=scf.Al2
# echo " 2" > $base.xyz
# echo "#cell 8.0 8.0 4.0 p p p" >> $base.xyz
# echo "Al   0 0 -1." >> $base.xyz
# echo "Al   0 0  1." >> $base.xyz

# base=scf.Al-fcc
# echo " 4" > $base.xyz
# echo "#cell 4.0 4.0 4.0 p p p" >> $base.xyz
# echo "Al   -1.0 -1.0 -1.0" >> $base.xyz
# echo "Al    1.0  1.0 -1.0" >> $base.xyz
# echo "Al    1.0 -1.0  1.0" >> $base.xyz
# echo "Al   -1.0  1.0  1.0" >> $base.xyz

# base=scf.Ne-fcc
# echo " 1" > $base.xyz  
# echo "%fcc_unitcell 0 2.235 2.235  2.235 0 2.235  2.235 2.235 0" >> $base.xyz
# echo "Ne  0 0 0" >> $base.xyz


### Cu LDA lattice constant from PHYSICAL REVIEW B 79, 085104 􏰀(2009􏰁), al. et Blaha
# base=scf.Cu-fcc
# echo " 4" > $base.xyz
# echo "#cell 3.522 3.522 3.522 p p p" >> $base.xyz
# echo "Cu   -.8805 -.8805 -.8805" >> $base.xyz
# echo "Cu    .8805  .8805 -.8805" >> $base.xyz
# echo "Cu    .8805 -.8805  .8805" >> $base.xyz
# echo "Cu   -.8805  .8805  .8805" >> $base.xyz

### Au LDA lattice constant 4.065 from Wikipedia
# base=scf.Au-fcc
# echo " 4" > $base.xyz
# echo "#cell 4.064 4.064 4.064 p p p" >> $base.xyz
# echo "Au   -1.016 -1.016 -1.016" >> $base.xyz
# echo "Au    1.016  1.016 -1.016" >> $base.xyz
# echo "Au    1.016 -1.016  1.016" >> $base.xyz
# echo "Au   -1.016  1.016  1.016" >> $base.xyz

# ### diamond LDA lattice constant 3.536 Ang from PHYSICAL REVIEW B 79, 085104 (2009), al. et Blaha
# base=scf.C-diamond
# echo " 8" > $base.xyz
# echo "#cell 3.536 3.536 3.536 p p p" >> $base.xyz
# # 3.536 / 8 == 0.442, 0.442 * 3 == 1.326
# echo "C  -1.326 -1.326 -1.326" >> $base.xyz
# echo "C   0.442  0.442 -1.326" >> $base.xyz
# echo "C  -0.442 -0.442 -0.442" >> $base.xyz
# echo "C   1.326  1.326 -0.442" >> $base.xyz
# echo "C  -1.326  0.442  0.442" >> $base.xyz
# echo "C   0.442 -1.326  0.442" >> $base.xyz
# echo "C  -0.442  1.326  1.326" >> $base.xyz
# echo "C   1.326 -0.442  1.326" >> $base.xyz

# base=scf.C_chain.sho
# echo " 1" > $base.xyz
# echo "#cell 1.420282 10 10 p p p" >> $base.xyz
# echo "C   0 0 0" >> $base.xyz

### silicon LDA lattice constant 5.378 Ang in diamond structure
# base=scf.silicon
# echo " 8" > $base.xyz
# echo "#cell 5.378 5.378 5.378 p p p" >> $base.xyz
### 5.378 / 8 == 0.67225, 0.67225 * 3 == 2.01675
# echo "Si  -2.01675 -2.01675 -2.01675" >> $base.xyz
# echo "Si   0.67225  0.67225 -2.01675" >> $base.xyz
# echo "Si  -0.67225 -0.67225 -0.67225" >> $base.xyz
# echo "Si   2.01675  2.01675 -0.67225" >> $base.xyz
# echo "Si  -2.01675  0.67225  0.67225" >> $base.xyz
# echo "Si   0.67225 -2.01675  0.67225" >> $base.xyz
# echo "Si  -0.67225  2.01675  2.01675" >> $base.xyz
# echo "Si   2.01675 -0.67225  2.01675" >> $base.xyz


### generate a control file
cat > control.sh << EOF
################################################################
#
#   a43 control file
#
################################################################
# use # for silent comments and #! for comments in the log file

## general verbosity of the log file
# verbosity=7

## display energies in custom units {Ha, Ry, eV}
#output.energy.unit=eV
## display distances in custom units {Bohr, Ang, nm, pm}
#output.length.unit=Ang

## atomic geometry
geometry.file=$base.xyz

## grid spacing or number of grid points of the dense grid
# basis=grid
grid.spacing.unit=Ang

## max. number of self-consistency iterations
self_consistency.max.scf=1
self_consistency.mix.density=0.25
atomic.valence.decay=0

## compute the Fermi level {exact, linearized}
# fermi.level=exact

## analyze the potentials up to vtot (DEBUG)
# potential_generator.use.bessel.projection=0
# potential_generator.use.direct.projection=0

## Poisson equation solver {fft, multigrid, none, cg, sd}
# electrostatic.solver=fft
# electrostatic.compensator=factorizable
electrostatic.compensator=generalized_Gaussian


##
## configuration of atomic PAW setups
##
## choose local potential method {V=parabola, V=sinc}

element_H="1s* 1 0 2p | 0.9 sigma .25 V=parabola"
element_C="2s* 2 2p* 2 0 3d 1e-100 | 1.2 sigma .43 V=parabola"

element_Al="3s* 2 3p* 1 0 3d | 2.05 sigma .645"
element_Si="3s* 2 3p* 2 0 3d | 2.0 sigma .643"
element_P=" 3s* 2 3p* 3 0 3d | 1.8 sigma .512"

## partial wave method {recreate_second, energy_ordering, classical, ...}
# single_atom.partial.wave.method=recreate_second

## relax partial wave in every SCF iteration {1:yes, 0:no}
single_atom.relax.partial.waves=0

## special verbosity for PAW setup
single_atom.init.echo=0

## special verbosity for PAW update
single_atom.echo=0

## bit mask for the first 50 atoms, -1:all, 1:only atom#0, 5:atoms#0 and #2 but not #1, ...
single_atom.echo.mask=3

## limit the number of partial waves per ell-channel, default=2
# single_atom.nn.limit=2

## optimize sigma for the occupied projectors {0:no, 1:yes, -1:optimize and show but don't use it, -10:optimize and exit}
# single_atom.optimize.sigma=1

## define via the localization of a state if it should be in the core
single_atom.core.state.localization=0.125

## debug options
# single_atom.synthetic.density.matrix=1
# single_atom.init.max.scf=0

## export PAW data in paw_xml format for GPAW or ABINIT
# single_atom.export.xml=0

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
bands.per.atom=4

## sampling of the Brillouin zone
# hamiltonian.kmesh.echo=9
# hamiltonian.kmesh=0
hamiltonian.kmesh=1
# hamiltonian.kmesh.x=21
## hamiltonian.floating.point.bits=64
hamiltonian.floating.point.bits=32

## configuration for basis=grid
# method of the grid eigensolver {cg, Davidson, none, explicit}
# grid.eigensolver=none
# grid.eigensolver=explicit
# grid.eigensolver=cg
conjugate_gradients.max.iter=4
grid.eigensolver.repeat=9

## for start wave functions use SHO functions with larger sigma spread
start.waves.scale.sigma=6

## load start waves from file (basis=grid), store wave functions to file
# start.waves=$base.waves.dat
store.waves=$base.waves.dat

## configuration for basis=pw
# plane_wave.solver {auto, both, direct, iterative}
# plane_wave.solver=direct
# plane_wave.cutoff.energy=520  eV
# plane_wave.cutoff.energy.unit=eV
# plane_wave.cutoff.energy=5.78  Ha

# plane_wave.solver=iterative
# plane_wave.iterative.solver.ratio=8.0
# plane_wave.iterative.solver=cg
# conjugate_gradients.max.iter=2
# plane_wave.max.cg.iterations=12
# davidson_solver.max.iterations=7

## also compute the eigenvalues of the overlap matrix?
# dense_solver.test.overlap.eigvals=0

# Export the Hamiltonian {0:not, 1:yes, -1:yes+abort} in format {xml, json}
hamiltonian.export=1
# hamiltonian.export.format=json
hamiltonian.export.format=xml

# Make structure_solver produce the same as potential_generator
# structure_solver.complex=1

# show variables defined in control (1=minimal, 2=unused, 4=defaults, 6=4+2, negative for details)
control.show=-1

EOF
## end of control file


for spacing in `seq 2 1 0`; do
  ## real-space grid spacing in Bohr
  project=$base.grid$spacing
  echo "# start calculation $project" && \
  $exe -test self_consistency \
        +control.file=control.sh \
        +basis=grid \
        +grid.spacing=`echo 0.222 / $spacing | bc -l` \
        "$@" > $project.out
  ./spectrum.sh $project.out > $project.spectrum.dat
done # spacing

for numax in `seq 9 2 9`; do
  ## SHO basis size numax
  project=$base.sho$numax
  echo "# start calculation $project" && \
  $exe -test self_consistency \
        +control.file=control.sh \
        +basis=sho \
        +grid.spacing=0.12501 \
        +sho_hamiltonian.test.numax=$numax \
        +sho_hamiltonian.use.sho_basis=yes \
        "$@" > $project.out
        # +sho_hamiltonian.test.sigma=1.0 \
  ./spectrum.sh $project.out > $project.spectrum.dat
done # numax

for ecut in `seq 11 3 1`; do
  ## plane-wave cutoff energy in Hartree
  project=$base.pw$ecut
  echo "# start calculation $project" && \
  $exe -test self_consistency \
        +control.file=control.sh \
        +basis=plane_wave \
        +plane_wave.cutoff.energy=$ecut \
        "$@" > $project.out
  ./spectrum.sh $project.out > $project.spectrum.dat
done # ecut
