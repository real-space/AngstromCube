**Data folder for AngstromCube**

[![DOI](https://img.shields.io/badge/GitHub-MIT-informational)](https://github.com/real-space/AngstromCube)

    AngstromCube uses both, the spherical representation and the Cartesian representation of
    Spherical Harmonic Oscillator (SHO) basis functions. To transform between the two representations,
    we can construct a unitary matrix from data in the *sho_unitary.dat* file.
    For details see Paul F. Baumeister and Shigeru Tsukamoto,
      [proceedings of PASC19](https://dl.acm.org/doi/10.1145/3324989.3325717)

**Files**
The root folder of this repository contains the following directories:
| File name       | Purpose                                                                        |
|-----------------|--------------------------------------------------------------------------------|
| sho_unitary.dat | transformation coefficients for 3D Spherical Harmonic Oscillator functions     |
| cho_unitary.dat | transformation coefficients for 2D Circular Harmonic Oscillator functions      |
