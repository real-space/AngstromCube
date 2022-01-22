#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

dir=~/Downloads/gpaw-setups-0.9.20000
ln -sf $dir ./gpaws

## see if the verbosity is monotonic
for f in `ls gpaws/*.LDA gpaws/*.PBE`; do

  echo -n "$f "
  $exe -t pawxml_import +pawxml_import.test.filename=$f

done

## see if we get the same after two conversions
f=C.LDA
$exe -t pawxml_import +pawxml_import.test.repeat=1 +pawxml_import.test.filename=$f
$exe -t pawxml_import +pawxml_import.test.repeat=1 +pawxml_import.test.filename=$f.repeat.xml
diff $f.repeat.xml $f.repeat.xml.repeat.xml
