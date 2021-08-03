#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43


#occ=auto
occ=custom

for Z in {1..120}; do

  out_Z=data.atom_core/atom_core.$Z.out
  echo -n "# " > $out_Z
  date >> $out_Z

# echo "# Z = $Z" ## show progress in terminal
  $exe +verbosity=5 \
    -test atom_core \
         +atom_core.test.Z=$Z \
         +atom_core.occupations=$occ \
       >> $out_Z

  grep ', E_tot= ' $out_Z
done
