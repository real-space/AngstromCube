#!/usr/bin/env bash

exe=./a43
project_base=green_dyadic

rm -rf $project_base.out
touch  $project_base.out

for nb in {1..20}; do
  echo -n "nb=$nb"
  $exe -test $project_base \
        +verbosity=7 \
        +$project_base.test.nb=$nb \
        "$@" > $project_base.out.1
  grep 'error' $project_base.out.1
  cat $project_base.out.1 >> $project_base.out
done
grep 'error' $project_base.out
