#!/usr/bin/env bash

exe=../src/a43

out=sho_potential.onsite.out
rm -f $out; echo -n "# " > $out; date >> $out

# for a in 0   1000    100 010 001     200 110 020 101 011 002      300 210 120 030 201 111 021 102 012 003   400 310 220 130 040 301 211 121 031 202 112 022 103 013 004; do
# for a in 0  1000    100 010 001     200 110 020 101 011 002  ; do
# for a in  200 020 002 ; do  ## the lowest non-linear ones
# for a in  200 ; do  ## the lowest non-linear one
# for a in 0 ; do  ## use a real local potential (not an artifical one)
# for a in   1000    100 010 001  ; do
# for a in 0 ; do
a=0; for p in {3..10}; do
# for a in 0 1000 100 200 300 400 ; do

(cd ../src/ && make -j) && \
  $exe +verbosity=9 \
    -test sho_potential. \
    +sho_potential.test.numax=2 \
    +sho_potential.test.sigma=2.0 \
    +sho_potential.test.method=6 \
    +sho_potential.test.sigma.asymmetry=1.1 \
    +sho_potential.test.artificial.potential=$a \
    +sho_potential.test.method4.percentage=$p \
    +sho_potential.test.show.potential.only=0 \
   >> $out

#     +sho_potential.test.lmax=8

done
grep '# V largest abs deviation of ' $out
