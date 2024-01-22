#!/usr/bin/env bash

exe=./a43

### name of the output file
out=pseudo_basis.xml

### write XML header
echo '<?xml version="1.0"?>'  > $out
echo '<basis version="0.1">' >> $out

### for all regular species
for Z in `seq 0 120`; do

    ### generate a pseudo_basis.$Z.xml file for each species
    $exe -t single_atom +single_atom.fit.basis.sigma=1.5 +single_atom.fit.basis.numax=19 +single_atom.fit.basis.weight=.0 +single_atom.test.Z=$Z

    ### remove the XML header and XML tail and append to pseudo_basis.xml
    sed -e '/<basis/d' -e '/<?xml/d' -e '/<\/basis>/d' pseudo_basis.$Z.xml >> $out

    ### delete the single files
    rm -f pseudo_basis.$Z.xml

done ### Z

### write XML tail
echo '</basis>' >> $out 
