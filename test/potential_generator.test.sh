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
# printf " 2 \n#cell 8 8 8 i i i \n" > $geometry_file
# echo "C  0 0 -0.65" >> $geometry_file
# echo "C  0 0  0.65" >> $geometry_file
### QE-spectrum          -19.6629 -10.5856  -7.2179  -7.2179  -7.1693  -0.4972  -0.4972  -0.1544 eV (self-consistent result)
# pg.C-dimer.sho3.out:#   -2.44596 12.4182 12.4182 15.8118 15.8118 17.2857 17.8086 26.7288 ... 227.187 274.533 eV
# pg.C-dimer.sho4.out:#   -10.4397 6.19268 9.4426 12.0828 12.0828 15.6831 15.6831 23.6335 ... 345.598 410.391 eV
# pg.C-dimer.sho5.out:#   -13.3748 3.22863 3.22863 4.31542 5.31541 8.59983 8.59983 20.4749 ... 473.253 556.877 eV
# pg.C-dimer.sho6.out:#   -16.5775 -0.639025 1.81662 2.97965 2.97965 8.39526 8.39526 18.6816 ... 617.558 713.264 eV
# pg.C-dimer.sho7.out:#   -17.4432 -1.76683 -1.44609 -1.44609 -0.128789 4.96216 4.96216 17.4456 ... 774.085 878.783 eV
# pg.C-dimer.sho8.out:#   -19.0737 -4.35561 -2.00145 -1.61998 -1.61998 4.75747 4.75747 16.2042 ... 940.895 1050.93 eV
# pg.C-dimer.sho9.out:#   -19.4077 -5.06115 -4.12805 -4.12805 -3.08344 2.8209 2.8209 15.5628 ... 1116.53 1228.39 eV
# pg.C-dimer.sho10.out#   -20.2971 -6.50843 -4.24552 -4.24552 -4.17229 2.64382 2.64382 14.6912 ... 1296.84 1411.54 eV
### exe using PW5 double -23.0606 -13.7356 -7.58652 -6.84787 -5.15872 -0.335201 -0.208443 0.065608 ... 133.323 133.324 eV
### exe using PW5  float -23.0604 -13.7355 -7.58635 -6.84766 -5.15843 -0.329002 -0.208301 0.0676805 ... 133.335 133.338 eV
### exe using PW10 float -23.6921 -14.7412 -8.86353 -8.05424 -5.99189 -1.53464 -1.03378 -0.233875 ... 269.28 269.28 eV
### exe using GRID float -26.633 -14.274 -13.599 -12.890 -10.141 -2.797 -2.786 -0.298 1.572 1.721 1.729 1.729 1.898 2.619 3.362 4.180 4.284 4.306 4.313 6.039 (did not converge)
# pg.C-dimer.pw1.out:#   -17.1491 -6.57399 -2.18672 -2.01646 -0.71439 1.13102 1.75672 2.21652 ... 25.2103 25.2193 eV
# pg.C-dimer.pw2.out:#   -20.5932 -10.4934 -4.36693 -3.86237 -2.20161 0.25182 1.69216 1.74257 ... 51.1207 51.1262 eV
# pg.C-dimer.pw3.out:#   -21.9315 -12.1795 -5.83802 -5.14459 -3.92636 -0.102485 1.10361 1.20521 ... 79.3067 79.315 eV
# pg.C-dimer.pw4.out:#   -22.6976 -13.2116 -6.96246 -6.24793 -4.78661 -0.225939 0.364226 0.441317 ... 107.346 107.358 eV
# pg.C-dimer.pw5.out:#   -23.0604 -13.7355 -7.58633 -6.8477 -5.15847 -0.328982 -0.208283 0.0677006 ... 133.335 133.338 eV
# pg.C-dimer.pw6.out:#   -23.3288 -14.1425 -8.11158 -7.31981 -5.46443 -0.714231 -0.539958 -0.120838 ... 161.676 161.677 eV
# pg.C-dimer.pw7.out:#   -23.4733 -14.3757 -8.41195 -7.59679 -5.64906 -1.03637 -0.729608 -0.1758 ... 189.534 189.539 eV
# pg.C-dimer.pw8.out:#   -23.5805 -14.5505 -8.64041 -7.8141 -5.81097 -1.26116 -0.866573 -0.20578 ... 213.086 213.086 eV
# pg.C-dimer.pw9.out:#   -23.6465 -14.6619 -8.78085 -7.9558 -5.91819 -1.42062 -0.963924 -0.222795 ... 243.497 243.497 eV
#            pw10        -23.6917 -14.7407 -8.86323 -8.05397 -5.99167 -1.53421 -1.03365 -0.23371 ... 269.281 269.281 eV

# project_base=pg.H-sc
# printf " 1 \n#cell 2.5 2.5 2.5 p p p \n" > $geometry_file
# echo "H  0 0 0" >> $geometry_file

# project_base=potential_generator.AlP
# printf " 2 \n#cell 4.233418 4.233418 8.466836 p p p \n" > $geometry_file
# printf " 2 \n#cell 10.5835 10.5835 12.7003 p p p \n" > $geometry_file
# printf " 2 \n#cell 21.16708996 21.16708996 25.400507952 p p p \n" > $geometry_file
# echo "Al   0 0 -1.058354498" >> $geometry_file
# echo "P    0 0  1.058354498" >> $geometry_file

# project_base=pg.Al-fcc
# printf " 4 \n#cell 4.0 4.0 4.0 p p p \n" > $geometry_file
# echo "Al   -1.0 -1.0 -1.0" >> $geometry_file
# echo "Al    1.0  1.0 -1.0" >> $geometry_file
# echo "Al    1.0 -1.0  1.0" >> $geometry_file
# echo "Al   -1.0  1.0  1.0" >> $geometry_file

# project_base=potential_generator.P-sc
# printf " 1 \n#cell 3.0 3.0 3.0 p p p \n" > $geometry_file
# echo "P 0 0 0" >> $geometry_file

# project_base=potential_generator.Al-sc
# printf " 1 \n#cell 3.0 3.0 3.0 p p p \n" > $geometry_file
# echo "Al 0 0 0" >> $geometry_file

# project_base=potential_generator.C_chain.sho
# printf " 1 \n#cell 1.420282 10 10 p p p \n" > $geometry_file
# echo "C   0 0 0" >> $geometry_file

# project_base=pg.He-atom
# printf " 1 \n#cell 6 6 6 p p p \n" > $geometry_file
# echo "He  0 0 0" >> $geometry_file

# project_base=pg.H-atom
# printf " 1 \n#cell 6 6 6 p p p \n" > $geometry_file
# echo "H  0 0 0" >> $geometry_file

# project_base=pg.Li-atom
# printf " 1 \n#cell 6 6 6 p p p \n" > $geometry_file
# echo "Li  0 0 0" >> $geometry_file


### generate a control file
cat > control.sh << EOF

# verbosity of the log files
verbosity=7

# show energies in units of electronVolt
output.energy.unit=eV
output.length.unit=Ang

# grid spacing of the dense grid (in Bohr)
potential_generator.grid.spacing=0.23622
#potential_generator.grid.spacing=0.1772   ## dense grid

# max number of self-consistency iterations
potential_generator.max.scf=1

# Poisson solver
electrostatic.solver=fft


# configuration of atomic PAW setups
#element_H="1s 1 0 | 0.9 sigma .308 V=parabola"
#element_He="1s 2 2p | 1.5 numax 1 sigma .537535"
#element_Li="2s 1 0 2p 2e-99 | 2.0 numax 1 sigma .7752 V=parabola"
#element_Li="2s 1 0 2p 2e-99 | 2.0 numax 2 sigma .612475 V=sinc"
#element_Li="2s 1 0 2p 2e-99 | 2.0 numax 1 sigma .8088 V=sinc"
#element_C="2s 2 2p 2 0 | 1.2 sigma .38 numax 2 V=sinc"
element_C="2s 2 2p 2 0 | 1.2 numax 1 sigma .445 V=parabola"
#element_C="2s 2 2p 2 0 | 1.2 numax 2 sigma .38 V=sinc"
#element_C="2s* 2 2p 2 0 3d | 1.2 numax 2 sigma .38 V=sinc"
#single_atom.partial.wave.energy=1.1 good for carbon with 2s*
#element_C="2s** 2 2p* 2 0 3d 4f | 1.2 numax 4 sigma .314327 V=sinc"
#element_C="2s*** 2 2p** 2 0 3d** 4f* 5g* 6h 7i | 1.2 numax 6 sigma .33055552 V=sinc"
#element_C="2s 2 2p 2 0 3d 4f 5g 6h 7i | 1.2 numax 6 sigma .33055552 V=sinc"
#element_C="2s 2 2p 2 0 | 1.2 numax 6 sigma .33055552 V=sinc"
#element_C="2s 2 2p 2 0 | 1.2 sigma .445 numax 9 V=parabola"
#element_Mg="2s 2 3s 2 2p 6 3p 2e-99 3d | 1.96 numax 3 sigma .60636 V=sinc"
#element_Mg="2s 2 3s* 2 2p 6 3p 2e-99 3d | 1.96 numax 4 sigma .578 V=sinc"
#element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma .5 V=parabola"
#element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1 V=sinc"
#single_atom.local.potential.method=sinc
single_atom.nn.limit=2
single_atom.partial.wave.method=energy_ordering
single_atom.echo=7
single_atom.init.echo=7
single_atom.optimize.sigma=0
single_atom.init.scf.maxit=0
# logarithmic derivatives
logder.unit=Ha
logder.start=-2
logder.stop=1
logder.step=1e-2

# configuration for basis=grid
bands.per.atom=10
# DEVEL option
devel.occupied.bands=3.5
eigensolver=cg
repeat.eigensolver=15
conjugate_gradients.max.iter=19
# for start wave functions use SHO functions with larger sigma spread
start.waves.scale.sigma=5
atomic.valence.decay=0
export.waves=waves.dat
#start.waves=waves.dat

# configuration for basis=sho
# spread of the SHO basis functions in Bohr
#sho_hamiltonian.test.sigma=2.0
#sho_hamiltonian.test.sigma=1.0
#sho_hamiltonian.test.sigma=0.7752
sho_hamiltonian.test.sigma=.5
#sho_hamiltonian.test.sigma=.308

# configuration for basis=sho or basis=pw
hamiltonian.test.kpoints=1
# start.waves.scale.sigma=1
hamiltonian.floating.point.bits=64

# configuration for basis=pw
pw_hamiltonian.density=1        0=no 1=yes
# pw_hamiltonian.solver= { auto both direct iterative }
pw_hamiltonian.solver=direct
davidson_solver.max.iterations=19
pw_hamiltonian.iterative.solver.ratio=2.0

# also compute the eigenvalues of the overlap matrix?
dense_solver.test.overlap.eigvals=1

# analyze the potentials up to vtot
potential_generator.use.bessel.projection=0
potential_generator.direct.projection=0

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

for ecut in `seq 7 1 7`; do
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
