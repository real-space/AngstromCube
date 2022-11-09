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

## compute the Fermi level {exact, linearized}
fermi.level=exact

## relax partial wave in every SCF iteration {1:yes, 0:no}
single_atom.relax.partial.waves=0

## special verbosity for PAW setup
single_atom.init.echo=3

## special verbosity for PAW update
single_atom.echo=3

## default for local potential method {parabola, sinc}
# single_atom.local.potential.method=sinc

## bit mask for the first 50 atoms, -1:all, 1:only atom#0, 5:atoms#0 and #2 but not #1, ...
single_atom.echo.mask=-1

## optimize sigma for the occupied projectors {0:no, 1:yes, -1:optimize and show but don't use it, -10:optimize and exit}
# single_atom.optimize.sigma=1

single_atom.core.state.localization=0.125

## debug options
# single_atom.synthetic.density.matrix=1
single_atom.init.max.scf=0

## export PAW data in paw_xml format for GPAW or ABINIT
single_atom.export.xml=0

###! the same radial grid is used to represent true and smooth quantities
# smooth.radial.grid.from=0

## to check logarithmic derivatives ensure start < stop
logder.unit=eV
logder.start=90
logder.stop=50
logder.step=1e-2

## show all used control environment variables at the end
control.show=3

