#!/usr/bin/env bash

echo -n "# " ## comment out the output of make, e.g. "make: Nothing to be done for `all'."
n=15
make -j &&
./a43 -t multi_grid. \
        +multi_grid.test.scheme=V \
        +multi_grid.test.pre.jacobi=2 \
        +multi_grid.test.post.jacobi=2 \
        +multi_grid.test.levels=$n \
        +multi_grid.test.min.level=1 \
        +multi_grid.test.iterations=30 \
        +multi_grid.test.rhs=sine \
        +verbosity=8 \
        | grep " smoothen level=$n 2 post" \
        | awk '{ print $9 }' >> mg.out
done

# make -j
# for n in {3..20} 
# do
# echo >> mg.out
# echo "# Number of levels is $n" >> mg.out
# ./a43 -t multi_grid. \
#         +multi_grid.test.scheme=V \
#         +multi_grid.test.pre.jacobi=2 \
#         +multi_grid.test.post.jacobi=2 \
#         +multi_grid.test.levels=$n \
#         +multi_grid.test.min.level=1 \
#         +multi_grid.test.iterations=200 \
#         +multi_grid.test.rhs=sine \
#         +verbosity=8 \
#         | grep " smoothen level=$n 2 post" \
#         | awk '{ print $9 }' >> mg.out
# done
