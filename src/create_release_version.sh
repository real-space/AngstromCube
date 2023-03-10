#!/usr/bin/env bash

### from https://www.muppetlabs.com/~breadbox/software/cppp.html
partial_preprocessor=~/Codes/cppp/cppp

defs=""\
" -D _VIEW2D_HAS_PARENTHESIS"\
" -D _VIEW3D_HAS_PARENTHESIS"\
" -D _VIEW3D_HAS_PARENTHESIS_2ARGS"\
" -D _VIEW3D_HAS_INDEXING"\
" -D _VIEW4D_HAS_PARENTHESIS"\
" -D _VIEW4D_HAS_PARENTHESIS_3ARGS"\
" -D _VIEW4D_HAS_INDEXING"

### hint: when moving a flag from defs to undefs, make sure to change -D to -U

undefs=""\
" -U NO_UNIT_TESTS"\
" -U DEBUG"\
" -U FULL_DEBUG"\
" -U DEVEL"\
" -U OLD_SINGLE_ATOM_UPDATE_INTERFACE"\
" -U LARGE_ANGULAR_GRIDS"\
" -U NEVER"

echo "Generation of a release version,"
echo "  apply the partial preprocessor"
echo "  with the following definitions"
echo "      $defs"
echo "  and the following undefs"
echo "      $undefs"

target_dir=../release
mkdir -p $target_dir

### delete the line that activates -D DEVEL in the development Makefile
sed '/DEVEL/d' Makefile > $target_dir/Makefile
### maybe some additional modifications of the FEATUREFLAGS in the Makefile
### will be necessary

echo "/*" > tmp_license_header_start
echo "*/" > tmp_license_header_end
for xxFile in `ls *.*xx *.h`
do
  # echo "$partial_preprocessor $defs $undefs $xxFile $target_dir"
    $partial_preprocessor $defs $undefs $xxFile $target_dir
  # Here, we could add a license header
  cat tmp_license_header_start ../LICENSE tmp_license_header_end $target_dir/$xxFile > tmp
  mv tmp $target_dir/$xxFile
done
rm tmp_license_header_*
