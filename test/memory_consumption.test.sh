#!/usr/bin/env bash

exe=./a43

bc_iii="i i i"
bc_pii="p i i"
bc_ppi="p p i"
bc_ppp="p p p"

# for bcs in bc_iii bc_pii bc_ppi bc_ppp; do
bcs="i i i"
# bcs="p p p"

for R in 10; do
  echo
  echo "# BoundaryCondition= $bcs Truncation Radius= $R"
  for ng in {1..24}; do
    xyz=empty_cell.xyz
    echo "0" > $xyz
    echo "#cell $ng $ng $ng $bcs" >> $xyz

    out=memory_consumption.r$R.n$ng.out

    $exe  -t parallel_potential \
          +verbosity=7 \
          +green_memory.show=1 \
          +geometry.file=$xyz \
          +grid.points=$((8*$ng)) \
          +grid.spacing=1.0 \
          +green_function.truncation.radius=$R \
          +green_function.sources=1 \
          +green_function.potential.exchange=0 \
          +energy_contour.bottom=0 \
          +energy_contour.parallel=1 \
          +energy_contour.fermidirac=0 \
          +energy_contour.matsubara=0 \
          +hamiltonian.kmesh=1 \
          +green_solver.iterations=0 \
          +control.show=-6 \
          > $out

    echo -n "$ng   "
    grep 'rank#0 green_memory now=' $out
  done # ng

done # R

# done # bcs