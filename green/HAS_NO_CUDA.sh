#!/usr/bin/env bash

for cuFile in `ls *.cu`
do
  fileName=${cuFile%.cu}
  ## create a soft link
  ln -s $cuFile ./$fileName.cxx
done
