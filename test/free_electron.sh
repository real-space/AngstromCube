#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

verbosity=8
for structure in fcc bcc sc ; do
        echo "free-electron for $structure, compute reference"
        for numax in {9..9} ; do ## show reference free electron parabolas
            out=free_electron.$structure.ref.out
            $exe -test sho_overlap. \
                      +sho_overlap.select.test=128 \
                      +sho_overlap.lattice.constant=8.0 \
                      +sho_overlap.crystal.structure=$structure \
                      +sho_overlap.test.DoS=1 \
                      +sho_overlap.test.Ref=1 \
                      +sho_overlap.kmesh.sampling=32 \
                      +sho_overlap.crystal.numax=$numax \
                      +verbosity=$verbosity > $out
        done
        for numax in {0..9} ; do
            echo "free-electron for $structure, compute numax=$numax"
            out=free_electron.$structure.nu$numax.out
            $exe -test sho_overlap. \
                      +sho_overlap.select.test=128 \
                      +sho_overlap.lattice.constant=8.0 \
                      +sho_overlap.crystal.structure=$structure \
                      +sho_overlap.test.DoS=1 \
                      +sho_overlap.test.Ref=0 \
                      +sho_overlap.kmesh.sampling=32 \
                      +sho_overlap.crystal.numax=$numax \
                      +verbosity=$verbosity > $out
        done
done
