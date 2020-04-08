#!/usr/bin/env bash

echo -n "# " ## comment out the output of make, e.g. "make: Nothing to be done for `all'."
make -j && \
./a43 -t multi_grid. \
        +multi_grid.test.scheme=V \
        +multi_grid.test.pre.jacobi=2 \
        +multi_grid.test.post.jacobi=2 \
        +multi_grid.test.levels=7 \
        +multi_grid.test.min.level=1 \
        +multi_grid.test.iterations=25 \
        +multi_grid.test.rhs=sine \
        +verbosity=8 \
#         | grep ' smoothen level=7 2 post' \
#         | awk '{ print $9 }'

# exit
# 
# ./a43 -t multi_grid. \
#         +multi_grid.test.scheme=V \
#         +multi_grid.test.pre.jacobi=2 \
#         +multi_grid.test.post.jacobi=2 \
#         +multi_grid.test.levels=15 \
#         +multi_grid.test.min.level=1 \
#         +multi_grid.test.iterations=15 \
#         +multi_grid.test.rhs=random \
#         +verbosity=8
#          | grep ' level = 15 after 2 post' | awk '{ print $12 }'
#         +verbosity=8 | awk '{ print $12 }'   ### show the final residual
#       +verbosity=8 > multi_grid.out.V_cycle.txt
##      +verbosity=20 > multi_grid.out.V_cycle.details
