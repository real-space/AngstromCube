#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

geometries="AlP_dimer.xyz Al_spherical_atom.xyz Cu40Zr22.xyz Cu_fcc_inversion.xyz Fe_bcc.xyz Sc2O3.xyz dna.xyz new_GST_set2.xyz tpa_on_gold.xyz
Al_single.xyz Cu320Zr180.xyz Cu_fcc.xyz Cu_fcc_shifted.xyz Po_sc.xyz graphene.xyz silicene.xyz"

referencedir=../ref
targetdir=.

for geo in $geometries ; do

    # run an A43 geometry analysis
    $exe +verbosity=7 \
      -test geometry_analysis \
      +geometry.file=geo/$geo \
      "$@" \
      > $targetdir/geometry_analysis.out_$geo

    # remove the line that displays the git key since that one will differ
    sed -ie '/git checkout/d' $targetdir/geometry_analysis.out_$geo
    # remove these lines to compare if order-N and N^2-algorithm produce the same
    sed -ie '/M atom-atom pairs/d' $targetdir/geometry_analysis.out_$geo

    # compare the ouput file with its reference version
#   diff     $targetdir/geometry_analysis.out_$geo $referencedir/geometry_analysis.out_$geo > $targetdir/geometry_analysis.dif_$geo
    diff -sq $targetdir/geometry_analysis.out_$geo $referencedir/geometry_analysis.out_$geo

done
grep 'recorded warning' $targetdir/geometry_analysis.out_* -A1

# cleanup
rm -f $targetdir/geometry_analysis.out_*.xyze
# rm -f $targetdir/geometry_analysis.out_*
# rm -f $targetdir/geometry_analysis.dif_*
