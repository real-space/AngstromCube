This directory stores the radial potential Zeff.00*Z* files
from a self-consistent generation of effective atom potentials.

A Zeff.00*Z* file is named with *Z* according to the nuclear number,
i.e. the number of protons in the core. For example, the file for
copper is named Zeff.029

The file contains two columns:
  - left  column: the radial coordinate *r* in Bohr
  - right column: the value *-rV(r)*, i.e. the potential multiplied by *-r* to remove the singularity

The user may confirm that the value at *r*=0 Bohr always matches *Z* exactly.
