#!/usr/bin/env bash

make -j -C ../src/
exe=../src/a43

### Alternative: compile standalone
# exe=./a.out
# sed -e 's/MODULE/load_balancer/g' test_MODULE.cxx > load_balancer_main.cxx
# g++ -std=c++11 -O0 -g -pedantic -Wall -I../include ../src/control.o ../src/recorded_warnings.o ../src/load_balancer.cxx load_balancer_main.cxx -o $exe
# rm -f load_balancer_main.cxx

project_base=load_balancer

rm -f $project_base.out
touch $project_base.out

### basic test: does the load_balancer produce equal distributions for exactly divisible numbers?
for n in {1..17}; do
  for m in {1..17}; do
    nm=$(( $n * $m ))
    echo -n "nprocs=$n x $m"
    $exe -test $project_base \
          +verbosity=7 \
          +$project_base.test.nprocs=$nm \
          +$project_base.test.nx=$n \
          +$project_base.test.ny=$m \
          "$@" | grep 'per process \[1, 1.00 +/- 0.00, 1\]' | wc -l
  done
done

exit


for nprocs in {1..116}; do
# echo -n "nprocs=$nprocs"
  $exe -test $project_base \
        +verbosity=7 \
        +$project_base.test.nprocs=$nprocs \
        +$project_base.test.nx=117 \
        +$project_base.test.ny=29 \
        +$project_base.test.nz=1 \
        "$@" >> $project_base.out
done
