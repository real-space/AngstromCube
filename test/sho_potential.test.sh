#!/usr/bin/env bash

exe=../src/a43

out=sho_potential.out
rm -f $out; echo -n "# " > $out
date >> $out

# method=1 # 'numerical'
# method=2 # 'between'
# method=3 # compare between to 'numerical'
# method=4 # 'on_site'
method=5 # compare 'on_site' to 'numerical'
# method=6 # compare 'on_site' to 'between'
# method=7 # compare 'between' and 'on_site' to 'numerical'
# method=-1 # same as method=7

make -j -C ../src

for percentage in {0..111} ; do # percentage of kinetic energy in method 'onsite'

    # art=1000 # an artifical constant potential
    # art=000 # a real potential loaded from +sho_potential.test.vtot.filename=vtot.dat

    for art in 000 ; do
#     for art in 000  1000  100 010 001 \
#                           200 110 020 101 011 002 \
#                           300 210 120 030 201 111 021 102 012 003 \
#                           400 310 220 130 040 301 211 121 031 202 112 022 103 013 004 ; do
    # for art in 0  1000    100 010 001     200 110 020 101 011 002 ; do
    # for art in  200 020 002 ; do  ## the lowest non-linear ones
    # for art in  200 ; do  ## the lowest non-linear one
    # for art in 0 ; do  ## use a real local potential (not an artifical one)
    # for art in   1000    100 010 001  ; do
    # for art in 0 ; do
    # art=0; for percentage in {3..10} ; do
    # for art in 0 1000 100 200 300 400 ; do
    
        out_p="sho_potential.out.percentage"

#         echo -n "# "
#         (cd ../src/ && make -j) && \
          $exe +verbosity=11 \
            -test sho_potential \
            +sho_potential.test.numax=3 \
            +sho_potential.test.sigma=2.0 \
            +sho_potential.test.method=$method \
            +sho_potential.test.sigma.asymmetry=1.1 \
            +sho_potential.test.artificial.potential=$art \
            +sho_potential.test.method.on_site.percentage=$percentage \
            +sho_potential.test.show.potential.only=1 \
          > $out_p

        # V largest abs deviation of 'between' to 'numerical' is (diag) 5.72209e-13 and 3.86358e-13 (off-diag), pot=000
        # V largest abs deviation of 'onsite' to 'numerical' is (diag) 1.51044e-10 and 0.0192199 (off-diag), pot=000
      
    done # art
    echo -n "$percentage % "
    grep '# V largest abs deviation of ' $out_p

    cat $out_p >> $out

done # percentage
