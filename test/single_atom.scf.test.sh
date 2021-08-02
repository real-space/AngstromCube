#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

for Z in {6..6}; do


    for method_name in atom_core single_atom; do
    
        ## compute the self-consistent solution using atom_core
        ## and
        ## compute the self-consistent solution using single_atom

        result_file=$method_name.scf.$Z.out
        $exe +verbosity=8 \
            -t $method_name \
              +$method_name.test.Z=$Z \
              +single_atom.nn.limit=1 \
              +single_atom.test.atomic.valence.density=1 \
              +single_atom.relax.partial.waves=0 \
              +single_atom.init.max.scf=66 \
              +single_atom.init.scf.echo=1 \
              +element_C="2s* 1.5 2p* 2.5 0 3d | 1.2 sigma .43" \
              +element_Mg="3s* 1.5 3p .5 3d | 1.96 sigma .41" \
              +element_Cn="7s 1.5 7p .5 6d 10 | 2. sigma .5" \
              +control.show=1 \
              "$@" > $result_file
              
      done
          
done

## Example magnesium:
#
# +element_Mg="3s* 2 3p 2e-99 3d | 1.96 sigma .41"
# atom_core   -199.454315121 Ha
# single_atom -199.451741448
# diff           0.002573673
#
# +element_Mg="3s* 1.5 3p .5 3d | 1.96 sigma .41"
# atom_core   -199.389617141 Ha
# single_atom -199.387011834
# diff           0.002605307
#
# diff diff    -31.634e-06 Ha
