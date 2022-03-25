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

echo "# real4 nu=0..5"
grep '# real4 nu=0 ' $project_base.out
grep '# real4 nu=1 ' $project_base.out
grep '# real4 nu=2 ' $project_base.out
grep '# real4 nu=3 ' $project_base.out
grep '# real4 nu=4 ' $project_base.out
grep '# real4 nu=5 ' $project_base.out

echo "# real8 nu=0..5"
grep '# real8 nu=0 ' $project_base.out
grep '# real8 nu=1 ' $project_base.out
grep '# real8 nu=2 ' $project_base.out
grep '# real8 nu=3 ' $project_base.out
grep '# real8 nu=4 ' $project_base.out
grep '# real8 nu=5 ' $project_base.out
