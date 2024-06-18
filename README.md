**AngstromCube**

[![DOI](https://img.shields.io/badge/GitHub-MIT-informational)](https://github.com/real-space/tfQMRgpu)

    AngstromCube is an experimental all-electron DFT application.
    Using the Green functions formalism, near-sightedness allows
    for linear-scaling and the projector augmented wave method (PAW), 
    in particular the revised PAW method (see Paul F. Baumeister and Shigeru Tsukamoto,
      proceedings of PASC'19, https://dl.acm.org/doi/10.1145/3324989.3325717)

![A4cube logo](doc/fig/a43_logo_bold.png)

**Name**
The name refers to a cube with edge length 1 Angstrom
which is abbreviated \AA in TeX code.
This is because always 4x4x4 real-space grid points are grouped
for performance which corresponds to roughly one \AA^3

**Principles**
The idea is to have a code that 
- is highly parallel
- can make use of GPUs
- does not require more input than the atomic coordinates
- can scale linearly

**Current Status**
- These features are ready:
    - MPI parallelization of parallel_potential and tfQMRgpu Green function solver
    - total energy calculation
    - complex wave functions, complex Green functions, ***k***-points
    - boundary conditions for Wfs: periodic and isolated, for Gf also repeat and vacuum
    - non-magnetic potential generation
    - MPI parallel Poisson solver for the electrostatics (no preconditioner)
    - SHO-projector PAW with all-electron atoms (currently only non-magnetic)
- These features are planned but have so far not been addressed:
    - different versions of LDA, GGA, meta-GGA (currently only LDA implemented)
    - efficient eigensolver for the grid Hamiltonian (currently inefficient subspace rotation method)
    - OpenMP parallelization (currently none)
    - GPU acceleration (currently none)
    - forces (currently none)
    - self-consistency convergence criteria (currently we set the number of iterations)
    - magnetism, collinear and non-collinear (currently only non-magnetic)
- Some features are build in only for development purposes:
    - a stable FFT Poisson solver for the electrostatic problem (serial only)
    - plane wave basis set using a dense matrix eigensolver (LAPACK) or iterative (in development)
    - dense eigensolver for the real-space grid Hamiltonian (expensive)
- These features are not intended to be implemented ever:
    - strain calculation
    - exact exchange
    - phonons

**Directories**
The root folder of this repository contains the following directories:
| Directory    | Purpose                                                                        |
|--------------|--------------------------------------------------------------------------------|
| data         | matrix element files for SHO transforms between radial and Cartesian bases     |
| doc          | documentation folder including manual and theory notes                         |
| external     | put third party libraries here                                                 |
| include      | source folder for C and C++ header files                                       |
| interfaces   | examples for the libliveatom.so library in C, Fortran90, Julia and Python      |
| ref          | reference outputs of certain unit tests                                        |
| src          | source folder for C++ and CUDA C++ sources                                     |
| test         | test scripts for certain modules                                               |
| tools        | experimental scripts in Julia and Rust (programming languages)                 |

**Abbreviations**
| Abbr. | Explanation                                                                           |
|-------|---------------------------------------------------------------------------------------|
| DFT   | Density Functional Theory                                                             |
| XC    | Exchange-correlation                                                                  |
| LDA   | Local Density Approximation                                                           |
| GGA   | Generalized Gradient Approximation                                                    |
| PAW   | Projector Augmented Wave                                                              |
| CPU   | Central Processing Unit                                                               |
| GPU   | Graphical Processing Unit                                                             |
| SHO   | Spherical Harmonic Oscillator                                                         |
| MPI   | Message Passing Interface                                                             |
| FFT   | Fast Fourier Transform                                                                |
| OMP   | OpenMP, Open Multi-Processing                                                         |
| TeX   | typesetting                                                                           |
| Gf    | Green function                                                                        |
| Wf    | Wave function (eigenstate)                                                            |
