**A4cube**

    A4cube is an all-electron DFT application under construction.
    It is based on a various representations: real-space grids and factorizable local orbitals.
    Furthermore, the Green functions formalism, near-sightedness and
    PAW, in particular revPAW
    (see Paul F. Baumeister, Shigeru Tsukamoto, proceedings of PASC19)
    
![A4cube logo](https://gitlab.version.fz-juelich.de/pbaum/a43/blob/devel/doc/fig/a43_logo_bold.png)

**Name**
    The name refers to a cube with edge length 1 Angstrom
    which is abbreviated \AA in LaTeX code.
    This is because always 4x4x4 real-space grid points are grouped
    for performance which corresponds to roughly one \AA^3
    
**Principles**
    The idea is to have a code that is 
	- highly parallel
    	- can make use of GPUs
    	- does not require more input than the atomic coordinates

**Current Status**
    These features are planned but have so far not been addressed:
    - different versions of LDA, GGA, meta-GGA (currenly only one LSDA: PZ81)
    - efficient eigensolver for the grid Hamiltonian (currently inefficient subspace rotation method)
    - MPI parallelization
    - OpenMP parallelization
    - GPU acceleration (currently none)
    - complex wave functions, k-points
    - fancy boundary conditions
    - total energy
    - forces
    - self-consistency (so far on the radial grids only)
    - magnetism
    These features are ready:
    - potential generation (non-magnetic)
    - multigrid Poisson solver for the electrostatics (serial)
    - SHO-projector PAW with all-electron atoms (non-magnetic)
    Some features are build in only for development:
    - FFT Poisson solver for the electrostatic problem
    - plane wave basis set using a dense eigensolver
    These feautures are not intended to be implemented ever:
    - Strain
    - phonons

**Reuse existing code modules**
    In order not to re-program and validate everything, we
    try to re-use code modules from 
        juRS
            MPI communicators and MPItools
            generation of the potential on the grid and inside the spheres
            input/output of the grid potential
            reading of positions
            revPAW implementation mod_all_electron.F90 
        ReSpawN
            Lebedev-Laikov grids angular_grid.mod
            revPAW implementation hermite_projectors.mod
            all-electron atom solver atom_core.mod
            gpu_bench/ - finite differences on GPUs
            sho_transform/ - projection and addition on GPUs
            unitary_transform.dat (only the data)
            unit_system.mod (now called display_units)
            energy_contour.mod (originally from KKRnano)
        miniKKR
            tfQMRgpu - linear solver on GPUs
        external:
            libxc?

**Reference version**
    In order to be able to test the code without GPUs,
    we need a module that replaces C, C++ parts
    and uses e.g. MKL to invert the Hamiltonian
    which is then stored explicitly 
    --> does not scale to bigger problems
        but useful to test correctness

**Local orbital method**
    For a method that is cheaper than the RS-grid based
    Green function version, we should implement
    a local-orbital method.
    Although, numerical atomic orbitals have shown the be
    a great basis in terms of physics, AAcube will have
    all the machinery to use PAW, so we can use pseudo
    atomic orbitals. Here, we can suggest another 
    SHO basis (with high \numax and small \sigma).
    \numax nSHO
0  1      okay combine 32 atoms
1  4      okay combine 8  atoms
2  10     -->  16 1/2 warp 37.5% , combine 3 atoms --> 30 of 32 6.3%
3  20     -->  32  1 warp  37.5% , combine 3 atoms --> 60 of 64 6.3%
4  35     -->  64  2 warps 45.3% pathological case
5  56     -->  64  2 warps 12.5%
6  84     -->  96  3 warps 12.5%
7  120    --> 128  4 warps  6.3%
8  165    --> 192  6 warps 14.1%
9  220    --> 224  7 warps  1.8%
10 286    --> 288  9 warps  0.7%
11 364    --> 384 12 warps  5.2%
12 455    --> 480 15 warps  5.2%
13 560    --> 576 18 warps  2.8%
14 680    --> 704 22 warps  3.4%
15 816    --> 832 26 warps  1.9%
    For this, we need to construct Hamiltonian and overlap in this basis.
    SHO states --> Overlap and Kinetic energy are simple
    1D analytical evaluation with shifted polynomials (at least for comparison).
    Potential operator has two possible evaluations:
        - Make a SHOT of V(x,y,z) with high lmax,
          use the Spherical2Radial transform to get V_lm(r)
          use Gaunt coefficients to get V_lm_l'm'(r)
          apply them to the basis functions in radial representation
          transform the representation with Spherical2Radial and its transpose.
          create the overlap of V|ket> with <bra|
    Or much simpler:
        - Write a GPU kernel that runs over the non-local regions
          of each pair of atoms. There, we compute the Hermite1D on the fly
          and integrate H_nx(x)*H_ny(y)*H_nz(z) * V(x,y,z) * H_nx'(x)*H_ny'(y)*H_nz'(z)
          on the grid. This can be vectorized over nSHO.
          The H1D for a != a' belong to different centers, so we need to expand the H1D for both.
          Block sizes need to be investigated, maybe 8^3 blocks are good
          if the potential is originally evaluated on a 2x denser grid anyway.
    How to chose the \sigma for the basis functions?
    Classical basis set estimates for SHO (3D) is
        nSHO * 8/pi^2/R^3 (1-r^2/R^2)^(3/2)
    Where the length scale (R) is given by the classical return radius.
    This is determined by the point where the cutoff energy touches
    the potential parabola. With this formula, we can get an impression 
    of how inhomogeneous the basis set density is in all space.

    Furthermore, we need a method to evaluate the density on the real-space grid.
    After solving with tfQMRgpu (or ParILU ??), we find a Green function.
    Its imaginary part integrated over energies is interpreted as density matrix D.
    It is non-local and symmetric. Thus we need to evaluate for each grid point
      rho(x,y,z) = sum_ia_i'a' D_ia_i'a' * H_nx(x)*H_ny(y)*H_nz(z) * H_nx'(x)*H_ny'(y)*H_nz'(z)
    where only atoms a and a' contribute that are close enough to the position (x,y,z).
    The H1D for a != a' belong to different centers, so we need to expand the H1D for both.
    Again, we can work with bitmasks and create the multiplied bitmasks for (c,a,a') beforehand
    Here, c is the cube index. This can be vectorized over (x,y,z).

**Eigenvalue solving method**
    Maintaining a method supporting Green's function (linear system solver)
    and eigenfunctions (eigensolver) is complicated.
    In particular because the PAW setup needs to be constructed differently.
    However, it would be a nice playground for methodologial experiments.
    The probably most suitable eigensolver in grid representation is ChASE.
    In the local orbital method, that could be ELPA (?? or MAGMA for GPUs??)

