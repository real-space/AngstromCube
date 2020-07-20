#!/usr/bin/env bash

exe=../src/a43

out=sho_potential.out
rm -f $out; echo -n "# " > $out; date >> $out

for a in 0   1000    100 010 001     200 110 020 101 011 002      300 210 120 030 201 111 021 102 012 003   400 310 220 130 040 301 211 121 031 202 112 022 103 013 004; do
# for a in   0  1000    100 010 001     200 110 020 101 011 002  ; do
# for a in  200 020 002 ; do  ## the lowest non-linear ones
# for a in  200 ; do  ## the lowest non-linear one
# for a in 0 ; do  ## use a real local potential (not an artifical one)

## ToDo: check if sigma.asymmetry works
##       check deviations of method=5

(cd ../src/ && make -j) && \
  $exe +verbosity=11 \
    -test sho_potential. \
    +sho_potential.test.numax=3 \
    +sho_potential.test.sigma=2.0 \
    +sho_potential.test.method=3 \
    +sho_potential.test.sigma.asymmetry=1.5 \
    +sho_potential.test.artificial.potential=$a \
   >> sho_potential.out

    #grep '# V ai#1 aj#0 ' sho_potential.out
done
grep 'V largest abs deviation of between to numerical' $out
