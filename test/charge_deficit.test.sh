#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43


method=C
for Z in {1..120} ; do
  out_Z=data.charge_deficit/$method.z$Z.out

for energy in {1..40} ; do
  E=`echo $energy \* 0.1 | bc -l`
  out=$method.out
  ./a43 -t single_atom \
          +single_atom.nn.limit=2 \
          +single_atom.test.Z=$Z \
          +single_atom.partial.wave.energy=$E \
          +single_atom.partial.wave.method=$method \
          +single_atom.stop.after.charge.deficit.eigenvalues=1 \
          +verbosity=7 \
          "$@" > $out

# for method in c m e r 1 2 ; do
#   E=1
#   out=$method.out
#   ./a43 -t single_atom \
#           +single_atom.nn.limit=2 \
#           +single_atom.partial.wave.method=$method \
#           +logder.step=1e-3 \
#           +verbosity=7 \
#           "$@" > $out
#           +element_Cu="4s* 1 0 3d 10 | 2.0 sigma .61" \
#        +single_atom.partial.wave.energy=-8 \ ### use energy derivatives instead of excited states

  echo ""                                            >> $out_Z
# echo "# method= $method"                           >> $out_Z
  echo "# energy split= $E Hartree"                  >> $out_Z
  grep 'eigenvalues of the s-overlap operator' $out  >> $out_Z
  grep 'eigenvalues of the p-overlap operator' $out  >> $out_Z
  grep 'eigenvalues of the d-overlap operator' $out  >> $out_Z
  grep 'eigenvalues of the f-overlap operator' $out  >> $out_Z
  echo ""                                            >> $out_Z

done # energy

done # Z

### Idea: use classical partial waves and classical numerical projectors
###       and expand the latter in a sigma=0.5 high numax basis, then,
###       check the eigenvalues of the charge deficit operator.
