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
printf " 1 \n#cell 8 8 8 i i i \n" > $geometry_file
echo "C  0 0 0" >> $geometry_file

# project_base=pg.Cu-atom
# printf " 1 \n#cell 2 2 2 p p p \n" > $geometry_file
# echo "Cu  0 0 0" >> $geometry_file

# project_base=pg.Mg-atom
# printf " 1 \n#cell 6 6 6 p p p \n" > $geometry_file
# echo "Mg  0 0 0" >> $geometry_file

# project_base=pg.C-dimer
# printf " 2 \n#cell 8 8 8 p p p \n" > $geometry_file
# echo "C  -0.65 0 0" >> $geometry_file
# echo "C   0.65 0 0" >> $geometry_file
## test translational invariance
# echo "C  0 0 -0.525" >> $geometry_file
# echo "C  0 0  0.775" >> $geometry_file

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

### Cu LDA lattice constant from PHYSICAL REVIEW B 79, 085104 􏰀(2009􏰁), al. et Blaha
# project_base=pg.Cu-fcc
# printf " 4 \n#cell 3.522 3.522 3.522 p p p \n" > $geometry_file
# echo "Cu   -.8805 -.8805 -.8805" >> $geometry_file
# echo "Cu    .8805  .8805 -.8805" >> $geometry_file
# echo "Cu    .8805 -.8805  .8805" >> $geometry_file
# echo "Cu   -.8805  .8805  .8805" >> $geometry_file

### diamond LDA lattice constant 3.536 Ang from PHYSICAL REVIEW B 79, 085104 􏰀(2009􏰁), al. et Blaha

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
#element_Al="3s* 2 3p* 1 0 3d | 1.8 sigma .5 V=parabola"
#element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1 V=sinc"

### Copper: for the occupation setting 4s 1 3d 10, the d-states should be 1.5 eV below the s-state
### element_Cu="4s 2 4p 2e-99 3d 5 4 | 2.0 sigma .651 V=sinc" (seems to work well with GRID, but
###                          fails with PW, s-state below d-states, ss-density matrix entry tiny)
### element_Cu="4s 1 0 4p 2e-99 3d 10 | 2.2 numax 2 sigma .71857 V=sinc" --> s-state below d-states, as well



### generate a control file
cat > control.sh << EOF

# configuration of atomic PAW setups
#element_C="2s 2 2p 2 0 | 1.2 numax 1 sigma .43 V=parabola"
element_C="2s 2 2p 2 0 | 1.2 numax 2 sigma .38 V=sinc"
#element_Cu="4s 1 0 4p 2e-99 3d 10 | 2.2 numax 2 sigma .742455 V=parabola"

# verbosity of the log files
verbosity=7

# show energies in units of electronVolt
output.energy.unit=eV
output.length.unit=Bohr

# grid spacing of the dense grid (in Bohr)
potential_generator.grid.spacing=0.23622
#potential_generator.grid.spacing=0.208
#potential_generator.grid.spacing=0.1772   ## dense grid

# max number of self-consistency iterations
potential_generator.max.scf=1

# Poisson solver {mg, fft, none, cg, sd} and {MG, load, Bessel0} in development
electrostatic.solver=fft

#single_atom.local.potential.method=sinc
single_atom.nn.limit=2
single_atom.partial.wave.method=energy_ordering
single_atom.echo=0
single_atom.init.echo=0
### bit mask for the first 50 atoms, -1:all, 1:only atom#0, 5:atoms#0 and #2 but not #1, ...
single_atom.echo.mask=1
#single_atom.optimize.sigma=0
single_atom.init.scf.maxit=0

#smooth.radial.grid.from=0

# logarithmic derivatives
logder.unit=Ha
logder.start=2
logder.stop=1
logder.step=1e-4

bands.per.atom=10

# configuration for basis=grid
# method of the grid eigensolver {cg, Davidson, none}
grid.eigensolver=cg
grid.eigensolver.repeat=3
conjugate_gradients.max.iter=19
# for start wave functions use SHO functions with larger sigma spread
start.waves.scale.sigma=6
atomic.valence.decay=0
start.waves=waves.dat
store.waves=waves.dat

# configuration for basis=sho
# spread of the SHO basis functions in Bohr
sho_hamiltonian.test.sigma=.5

# configuration for basis=sho or basis=pw
hamiltonian.kmesh.x=1
# start.waves.scale.sigma=1
hamiltonian.floating.point.bits=64

# configuration for basis=pw
# pw_hamiltonian.solver {auto, both, direct, iterative}
pw_hamiltonian.solver=direct

#pw_hamiltonian.solver=iterative
davidson_solver.max.iterations=9
pw_hamiltonian.iterative.solver.ratio=4.0
pw_hamiltonian.iterative.solver=cg

# also compute the eigenvalues of the overlap matrix?
#dense_solver.test.overlap.eigvals=0

# analyze the potentials up to vtot
#potential_generator.use.bessel.projection=0
#potential_generator.direct.projection=0

single_atom.init.scf.maxit=1
#single_atom.export.xml=1

EOF


for spacing in `seq 1 1 1`; do
  project=$project_base.grid$spacing
  (cd ../src/ && make -j) && \
  echo "# start calculation $project" && \
  $exe -test potential_generator. \
        +control.file=control.sh \
        +basis=grid \
        $1 \
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
        $1 \
        > $project.out
        ./spectrum.sh $project.out > $project.spectrum.dat
done

for ecut in `seq 4 2 0`; do
  project=$project_base.pw$ecut
  (cd ../src/ && make -j) && \
  echo "# start calculation $project" && \
  $exe -test potential_generator. \
        +control.file=control.sh \
        +basis=pw \
        +pw_hamiltonian.cutoff.energy=$ecut \
        $1 \
        > $project.out
        ./spectrum.sh $project.out > $project.spectrum.dat
done

## results for C-dimer
### QE-spectrum          -19.6629 -10.5856  -7.2179  -7.2179  -7.1693  -0.4972  -0.4972  -0.1544 eV (self-consistent result)
# pg.C-dimer.sho3.out:#   -2.44596 12.4182 12.4182 15.8118 15.8118 17.2857 17.8086 26.7288 ... 227.187 274.533 eV
# pg.C-dimer.sho4.out:#   -10.4397 6.19268 9.4426 12.0828 12.0828 15.6831 15.6831 23.6335 ... 345.598 410.391 eV
# pg.C-dimer.sho5.out:#   -13.3748 3.22863 3.22863 4.31542 5.31541 8.59983 8.59983 20.4749 ... 473.253 556.877 eV
# pg.C-dimer.sho6.out:#   -16.5775 -0.639025 1.81662 2.97965 2.97965 8.39526 8.39526 18.6816 ... 617.558 713.264 eV
# pg.C-dimer.sho7.out:#   -17.4432 -1.76683 -1.44609 -1.44609 -0.128789 4.96216 4.96216 17.4456 ... 774.085 878.783 eV
# pg.C-dimer.sho8.out:#   -19.0737 -4.35561 -2.00145 -1.61998 -1.61998 4.75747 4.75747 16.2042 ... 940.895 1050.93 eV
# pg.C-dimer.sho9.out:#   -19.4077 -5.06115 -4.12805 -4.12805 -3.08344 2.8209 2.8209 15.5628 ... 1116.53 1228.39 eV
# pg.C-dimer.sho10.out#   -20.2971 -6.50843 -4.24552 -4.24552 -4.17229 2.64382 2.64382 14.6912 ... 1296.84 1411.54 eV
### exe using GRID       -21.702 -10.513 -9.352 -9.352 -8.465 -2.165 -2.165 -0.296 1.570 1.705 1.705 1.721 1.909 2.608 3.362 4.180 4.182 4.183 4.448 6.039
### exe using PW1..10 (double)
# 0.000000     -16.8874 -6.44077 -2.20068 -2.20068 -0.729291 1.11406 1.75765 2.24427    25.2062 25.2209
# 0.000000     -19.9523 -9.47263 -4.42102 -4.42102 -2.55295 0.190194 1.74456 1.77841    51.1208 51.1264
# 0.000000     -20.9086 -10.0109 -5.95291 -5.95291 -5.03678 -0.131462 1.11871 1.11871    79.3067 79.315
# 0.000000     -21.3335 -10.2347 -7.16013 -7.16013 -6.47815 -0.212975 0.161179 0.161179    107.345 107.347
# 0.000000     -21.5164 -10.3884 -7.85391 -7.85391 -7.0845 -0.509519 -0.509519 -0.238227    133.335 133.339
# 0.000000     -21.6531 -10.5269 -8.45097 -8.45097 -7.55846 -1.11633 -1.11633 -0.257763    161.676 161.677
# 0.000000     -21.7446 -10.6177 -8.79659 -8.79659 -7.85247 -1.50638 -1.50638 -0.269901    189.531 189.54
# 0.000000     -21.8238 -10.6678 -9.06692 -9.06692 -8.1441 -1.78376 -1.78376 -0.281995    213.086 213.087
# 0.000000     -21.8833 -10.7063 -9.23394 -9.23394 -8.34004 -1.9806 -1.9806 -0.290034    243.497 243.497
# 0.000000     -21.9348 -10.7368 -9.33212 -9.33212 -8.48694 -2.12323 -2.12323 -0.295576    269.281 269.281
