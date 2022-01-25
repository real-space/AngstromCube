#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43


result_file=single_atom.pawxml_import.out
rm -f $result_file
echo -n "# " > $result_file
date        >> $result_file

#for Z in {1..86}; do

### all existing Sy.LDA potentials
for Z in  1  2  3  4  5  6  7  8  9 10 \
         11 12 13 14 15 16 17 18 19 20 \
         21 22 23 24 25 26 27 28 29 30 \
         31 32 33 34 35 36 37 38 39 40 \
         41 42    44 45 46 47 48 49 50 \
         51 52 53 54 55 56 \
            72 73 74 75 76 77 78 79 80 \
         81 82 83       86
do

    ## experiment: try to load GPAW pawdata
    $exe +verbosity=7 \
    -test single_atom \
          +single_atom.test.Z=$Z \
          +single_atom.pawxml.path=gpaws \
          +single_atom.select.test=2 \
          >> $result_file

done # Z
