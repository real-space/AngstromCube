#!/usr/bin/env bash

partial_preprocessor=~/Codes/cppp/cppp

defs="-D NO_UNIT_TESTS"\
"-D _VIEW2D_HAS_PARENTHESIS -D _VIEW3D_HAS_PARENTHESIS"\
"-D _VIEW3D_HAS_PARENTHESIS_2ARGS -D _VIEW3D_HAS_INDEXING"\
"-D _VIEW4D_HAS_PARENTHESIS -D _VIEW4D_HAS_PARENTHESIS_3ARGS"\
"-D _VIEW4D_HAS_INDEXING"

undefs="-U DEBUG -U FULL_DEBUG -U DEVEL"\
"-U OLD_SINGLE_ATOM_UPDATE_INTERFACE"\
"-U LARGE_GRIDS"

echo "Generation of a release version,"
echo "  apply the partial preprocessor"
echo "  with the following definitions"
echo "      $defs"
echo "  and the following undefs"
echo "      $undefs"

source_dir=./
target_dir=$1

cp $source_dir/Makefile $target_dir
for xxFile in `ls $source_dir/*.*xx`
do
    $partial_preprocessor $defs $undefs $xxFile > $target_dir/$xxFile
done
