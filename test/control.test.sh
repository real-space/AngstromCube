#!/usr/bin/env bash

exe=./a43

out=control.test.out
rm -f $out; echo -n "# " > $out
date >> $out

for show in {-7..7} ; do

  echo -n "# "
  $exe -test control +verbosity=1 +control.show=$show >> $out

done

grep 'variables listed' $out
