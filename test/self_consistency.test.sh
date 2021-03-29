#!/usr/bin/env bash

exe=../src/a43
geometry_file=atoms.xyz

project_base=scf.C-atom
printf " 1 \n#cell 8 8 8 i i i \n" > $geometry_file
echo "C  0 0 0" >> $geometry_file


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
self_consistency.grid.spacing.unit=Ang
#self_consistency.grid.spacing=0.083334
#self_consistency.grid.spacing=0.10416667
self_consistency.grid.points=24

## max. number of self-consistency iterations
self_consistency.max.scf=1
self_consistency.mix.density=0.25
atomic.valence.decay=1
## compute the Fermi level {exact, linearized}
fermi.level=exact

## analyze the potentials up to vtot (DEBUG)
#self_consistency.use.bessel.projection=0
#self_consistency.use.direct.projection=0

## Poisson equation solver {multigrid, none, cg, sd, fft}
electrostatic.solver=fft


##
## configuration of atomic PAW setups
##

## configurations of elements
#element_C="2s* 2 2p 2 0 | 1.2 numax 2 sigma .4304 V=parabola"
element_C="2s* 2 2p* 2 0 | 1.2 numax 3 sigma .345353 V=parabola"

#element___="1s 1e-6 2p 0 | 1.0 numax 1 sigma .5 V=parabola"

## partial wave method {energy_ordering, recreate_second, classical, ...}
single_atom.partial.wave.method=recreate_second

## relax partial wave in every SCF iteration {1:yes, 0:no}
single_atom.relax.partial.waves=0

## special verbosity for PAW setup
single_atom.init.echo=0

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
single_atom.init.scf.maxit=0

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
hamiltonian.kmesh.echo=9
#hamiltonian.kmesh.x=3
hamiltonian.kmesh=0
hamiltonian.floating.point.bits=64

## configuration for basis=grid
# method of the grid eigensolver {cg, Davidson, none, explicit}
grid.eigensolver=none
#grid.eigensolver=explicit
#grid.eigensolver=cg
#conjugate_gradients.max.iter=1
#grid.eigensolver.repeat=1

## for start wave functions use SHO functions with larger sigma spread
start.waves.scale.sigma=6

## load start waves from file, store wave functions to file
#start.waves=$project_base.waves.dat
store.waves=$project_base.waves.dat

## configuration for basis=pw
#pw_hamiltonian.solver {auto, both, direct, iterative}
pw_hamiltonian.solver=direct

#pw_hamiltonian.solver=iterative
#pw_hamiltonian.iterative.solver.ratio=8.0
#pw_hamiltonian.iterative.solver=cg
#conjugate_gradients.max.iter=2
#pw_hamiltonian.max.cg.iterations=12
#davidson_solver.max.iterations=7

## also compute the eigenvalues of the overlap matrix?
#dense_solver.test.overlap.eigvals=0

export.hamiltonian=1            # 0:not, 1:yes, -1:yes+abort
# xml or json
export.hamiltonian.format=xml

EOF


for spacing in `seq 2 1 2`; do
  project=$project_base.grid$spacing
  (cd ../src/ && make -j) && \
  echo "# start calculation $project" && \
  $exe -test self_consistency \
        +control.file=control.sh \
        +basis=grid \
        $1 > $project.out
  ./spectrum.sh $project.out > $project.spectrum.dat
done
