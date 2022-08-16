#!/usr/bin/env bash

make -j -C ../src
exe=../src/a43



for Z in {6..6}; do
  for numax in `seq 1 5`; do
    result_file=single_atom.optimize_sigma.out
    rm -f $result_file
    touch $result_file
    for rcut in `seq -f "%.2f" 0.9 0.01 2.5`; do
        echo "# numax = $numax    rcut = $rcut Bohr    Z = $Z" ## show progress in terminal
        out_Z=single_atom.$Z.$numax.$rcut.out

        ## experiment: check the optimized sigma as a function of rcut, with parameter numax
        $exe +verbosity=3 \
        -test single_atom \
             +single_atom.test.Z=$Z \
             +single_atom.nn.limit=1 \
             +single_atom.optimize.sigma=-10 \
             +logder.start=2 +logder.stop=1 \
             +element_C="2s 2 2p 2 0 | $rcut numax $numax sigma .5 V=parabola" \
             > $out_Z

        echo -n "rcut= $rcut Bohr " >> $result_file
        grep 'optimized sigma=' $out_Z | grep "with quality" >> $result_file

    done
    awk '{print $2 " " $8 " " $16}' $result_file > sigma_of_rcut.z$Z.n$numax.dat
  done
done

# exit


echo "# lowest s-eigenvalues " > eigs_s0
echo "# higher s-eigenvalues " > eigs_s1
echo "# lowest p-eigenvalues " > eigs_p0
echo "# higher p-eigenvalues " > eigs_p1
echo "# lowest d-eigenvalues " > eigs_d0
echo "# lowest f-eigenvalues " > eigs_f0

### analyze
for Z in {1..47}; do
  echo
  echo "# Z = $Z" ## show progress in terminal
  out_Z=single_atom.$Z.out
  grep '#  lowest eigenvalues for ell=0 ' $out_Z
  grep "#  valence  [1-9]s " $out_Z
  grep '#  lowest eigenvalues for ell=1 ' $out_Z
  grep "#  valence  [1-9]p " $out_Z
  grep '#  lowest eigenvalues for ell=2 ' $out_Z
  grep "#  valence  [1-9]d " $out_Z

  echo -n "$Z " >> eigs_s0; grep '#  lowest eigenvalues for ell=0 ' $out_Z | awk '{print $6 }' >> eigs_s0
  echo -n "$Z " >> eigs_s1; grep '#  lowest eigenvalues for ell=0 ' $out_Z | awk '{print $7 }' >> eigs_s1
  echo -n "$Z " >> eigs_s0; grep "#  valence  [1-9]s " $out_Z | head -1    | awk '{print $6 }' >> eigs_s0
  echo -n "$Z " >> eigs_s1; grep "#  valence  [1-9]s " $out_Z | tail -1    | awk '{print $6 }' >> eigs_s1

  echo -n "$Z " >> eigs_p0; grep '#  lowest eigenvalues for ell=1 ' $out_Z | awk '{print $6 }' >> eigs_p0
  echo -n "$Z " >> eigs_p1; grep '#  lowest eigenvalues for ell=1 ' $out_Z | awk '{print $7 }' >> eigs_p1
  echo -n "$Z " >> eigs_p0; grep "#  valence  [2-9]p " $out_Z | head -1    | awk '{print $6 }' >> eigs_p0
  echo -n "$Z " >> eigs_p1; grep "#  valence  [2-9]p " $out_Z | tail -1    | awk '{print $6 }' >> eigs_p1

  echo -n "$Z " >> eigs_d0; grep '#  lowest eigenvalues for ell=2 ' $out_Z | awk '{print $6 }' >> eigs_d0
  echo -n "$Z " >> eigs_d0; grep "#  valence  [3-9]d " $out_Z | head -1    | awk '{print $6 }' >> eigs_d0

  echo -n "$Z " >> eigs_f0; grep '#  lowest eigenvalues for ell=3 ' $out_Z | awk '{print $6 }' >> eigs_f0
  echo -n "$Z " >> eigs_f0; grep "#  valence  [3-9]d " $out_Z | head -1    | awk '{print $6 }' >> eigs_f0
  
done
# exit



for Z in {1..47}; do
  ## a separate out file for each Z
  out_Z=single_atom.$Z.out
  echo -n "# using sigma_config    " > $out_Z
  date >> $out_Z

  echo "# Z = $Z" ## show progress in terminal
  $exe +verbosity=7 \
  -test single_atom \
       +single_atom.test.Z=$Z \
       +single_atom.test.atomic.valence.density=1 \
       +single_atom.nn.limit=2 \
       +single_atom.init.max.scf=1 \
       +single_atom.optimize.sigma=-1 \
       +element_H="1s* 1 0 2p 2e-99 | 0.9 sigma .623" \
       +element_He="1s* 2 2p | 1.5 sigma .75" \
       +element_B="2s* 2 2p* 1 0 3d | 1.2 sigma .7" \
       +element_C="2s* 2 2p* 2 0 3d | 1.2 sigma .8" \
       +element_N="2s* 2 2p* 3 0 3d | 1.0 sigma .8" \
       +element_O="2s* 2 2p* 3 1 3d | 1.13 sigma .9" \
       +element_F="2s* 2 2p* 3 2 3d | 1.2 sigma .8" \
       +element_Ne="2s* 2 2p* 6 3d | 1.8 sigma .7" \
       +element_Mg="2s 2 3s 2 2p 6 3p 2e-99 3d | 1.96 sigma .41" \
       +element_P="3s* 2 3p* 3 0 3d | 1.8 sigma 1.1" \
       +element_S="3s* 2 3p* 3 1 3d | 1.7 sigma 1." \
       +element_Ar="3s* 2 3p* 6 3d | 1.6 sigma .8" \
       +element_K="3s 2 4s 1 0 3p* 6 3d | 1.77 sigma .8" \
       +element_Sc="3s 2 4s 2 3p 6 4p 2e-99 3d* 1 0 | 2.32 sigma .6" \
       +element_Cu="4s* 1 1 4p* 2e-99 3d* 9 4f 2e-99 | 2.0 sigma .61" \
       +element_Ge="4s* 2 4p* 2 0 4d | 1.9 sigma 1." \
       +element_As="4s* 2 4p* 3 0 4d | 2.0 sigma 1.1" \
       +element_Se="4s* 2 4p* 3 1 4d | 1.6 sigma 1." \
       +element_Kr="4s* 2 4p* 6 4d | 2.2 sigma .9" \
       +element_Rb="4s 2 5s 1 0 4p* 6 4d | 2.3 sigma 1." \
       +element_u0="7s 2 8s 2 7p 6 8p 2e-99 6d 10 6f 2e-99 | 3. sigma .9" \
       +logder.start=2.5 +logder.step=.001 +logder.stop=1 \
      > single_atom.$Z.out

done

# exit








for Z in {1..120}; do
  ## a separate out file for each Z
  out_Z=single_atom.out.$Z.try
  echo -n "# using sigma_config    " > $out_Z
  date >> $out_Z

  echo "# Z = $Z" ## show progress in terminal
  $exe +verbosity=7 \
   -test single_atom \
        +single_atom.test.Z=$Z \
        +single_atom.config=custom \
        +single_atom.test.atomic.valence.density=1 \
        +single_atom.nn.limit=2 \
        +logder.start=-2 +logder.stop=-3 +logder.step=.001 \
        >> $out_Z
done
### when the energy difference =0 it seems like the eigenvalues of the charge deficit for ell=0 are always -1 and some positive value up to Z <= 86 and for Z=119,120

exit




### We find that when we use a fixed energy parameter, 
###     the eigenenergies found in eigenstate_analysis match almost exactly 
###     that energy parameter as long as be are bound (E < 0)

out=single_atom.out.conf
echo -n "# using sigma_config    " > $out
date     >> $out
echo " " >> $out

for Z in {117..117}; do
  ## either make a separate out file for each Z
  out_Z=single_atom.out.$Z.conf.new
  echo -n "# using sigma_config    " > $out_Z
  date >> $out_Z
  ## or concat all into the global out file
  # out_Z=$out

  echo "# Z = $Z" ## show progress in terminal
  $exe +verbosity=7 \
   -test single_atom \
        +single_atom.test.Z=$Z \
        +single_atom.partial.wave.method=m \
        +single_atom.partial.wave.energy.derivative=0 \
        +single_atom.partial.wave.energy=1.0 \
        +single_atom.test.atomic.valence.density=1 \
        +single_atom.from.sigma.config=1 \
        +single_atom.nn.limit=1 \
        +logder.start=-2 +logder.stop=1 +logder.step=.001 \
        +show.state.diagram=-1 \
        >> $out_Z
done
## We find stability execept for ell=2 of Z=7--9 and ell=0 of Z=80,84--92,117--120
##   which we need to analyze further for the "well-configured" elements (Z<=86)

