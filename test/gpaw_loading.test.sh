#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

dir=~/Downloads/gpaw-setups-0.9.20000
ln -sf $dir ./gpaws

## see if the verbosity is monotonic
for f in `ls gpaws/*.LDA gpaws/*.PBE`; do

  echo -n "$f "
  $exe -t gpaw_loading +gpaw_loading.filename=$f

done
