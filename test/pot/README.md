This directory stores the radial potential Zeff.00*Z* files
from a self-consistent generation of effective atom potentials.

A Zeff.00*Z* file is named with *Z* according to the nuclear number,
i.e. the number of protons in the core. For example, the file for
copper is named Zeff.029

The file contains two columns:
  - left  column: the radial coordinate *r* in Bohr
  - right column: the value *-rV(r)* in Bohr*Hartree. 
  
The potential *V(r)* is multiplied by *-r* to remove the singularity. 
Consequently, the values at *r=0* are equal to *Z* as the potential *V(r)= -Z/r* cannot be screened at the origin.
Without screening, all entries of the right column would read *Z*.
