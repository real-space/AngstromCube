#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

# $exe +verbosity=17 -test sigma_config \
#      +element_Lu="5s* 2 6s 2 5p 6 6p 2e-99 5d* 1 .1 | 2.4 sigma .6 Z 71.1 3pCore 3 2.5 " \
#      "$@" > sigma_config.out

     ### should launch this warning:
     ### sigma_config.cxx:549 warn("PAW setup for Lu (Z=71.1) is charged with -0.5 electrons")

dir=data.sigma_config
mkdir -p $dir

# for Z in {1..86}; do
for Z in {6..6}; do

    out=$dir/sigma_config.$Z.out
    rm -f $out
    touch $out

    # show the configuration string
    $exe -t sigma_config +sigma_config.test.Z=$Z "$@" | grep element_ >> $out
    
#   for m in m e C c r 2 1; do
    for m in m c r; do
        echo "#"                          >> $out
        echo "# Partial Wave method = $m" >> $out
    
        $exe  +verbosity=7 \
            -t single_atom \
              +single_atom.test.Z=$Z \
              +single_atom.optimize.sigma=-1 \
              +single_atom.partial.wave.method=$m \
              +logder.step=1e-2 \
              "$@" >> $out
#               +single_atom.export.xml=1 \
#               +single_atom.export.path=$dir \
#               +single_atom.smooth.radial.grid.from=0 \
        echo "# Partial Wave method = $m" >> $out
        echo "#"                          >> $out
        echo "#"                          >> $out
        echo "#"                          >> $out
        echo "#"                          >> $out
        echo "#"                          >> $out
        echo "#"                          >> $out
        echo "#"                          >> $out

    done

    # show the configuration string again
    $exe -t sigma_config +sigma_config.test.Z=$Z "$@" | grep element_ >> $out
    

done
