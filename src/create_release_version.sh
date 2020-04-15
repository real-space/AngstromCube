#!/usr/bin/env bash

### from https://www.muppetlabs.com/~breadbox/software/cppp.html
partial_preprocessor=~/Codes/cppp/cppp

defs=""\
" -D NO_UNIT_TESTS"\
" -D _VIEW2D_HAS_PARENTHESIS"\
" -D _VIEW3D_HAS_PARENTHESIS"\
" -D _VIEW3D_HAS_PARENTHESIS_2ARGS"\
" -D _VIEW3D_HAS_INDEXING"\
" -D _VIEW4D_HAS_PARENTHESIS"\
" -D _VIEW4D_HAS_PARENTHESIS_3ARGS"\
" -D _VIEW4D_HAS_INDEXING"

undefs=""\
" -U DEBUG"\
" -U FULL_DEBUG"\
" -U DEVEL"\
" -U OLD_SINGLE_ATOM_UPDATE_INTERFACE"\
" -U LARGE_GRIDS"

echo "Generation of a release version,"
echo "  apply the partial preprocessor"
echo "  with the following definitions"
echo "      $defs"
echo "  and the following undefs"
echo "      $undefs"

target_dir=../release

cp Makefile $target_dir
### maybe some modifications of the FEATUREFLAGS in the Makefile 
### will be necessary

for xxFile in `ls *.*xx *.h`
do
  # echo "$partial_preprocessor $defs $undefs $xxFile $target_dir"
    $partial_preprocessor $defs $undefs $xxFile $target_dir
done
