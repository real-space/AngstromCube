#!/usr/bin/env bash

exe=../src/a43
project_base=green_function

### generate a control file
cat > control.sh << EOF
################################################################
#
#   a43 control file
#
################################################################
# use # for silent comments and #! for comments in the log file

## general verbosity of the log file
verbosity=7

## display energies in custom units {Ha, Ry, eV}
output.energy.unit=eV
## display distances in custom units {Bohr, nm, Ang}
output.length.unit=Bohr

# green.function.truncation.radius=10
green.function.benchmark.iterations=-1

EOF

for rcut in `seq 0 0.1 9`; do
  project=$project_base.r$rcut
  (cd ../src/ && make -j) && \
#   echo "# start calculation $project" && \
  $exe -test green_function \
        +control.file=control.sh \
        +green.function.truncation.radius=$rcut \
        | grep 'target blocks per source block'
#         $1 > $project.out
done
