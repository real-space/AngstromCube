#!/usr/bin/env bash

grep 'spectrum  ' $1                             > tmp1
sed -e 's/spectrum/  /g'                    tmp1 > tmp2
sed -e 's/\.\.\./  /g'                      tmp2 > tmp3
sed -e 's/\# //g'                           tmp3 > tmp4
sed -e 's/ Ha//g' -e 's/ eV//g'             tmp4

rm -f tmp1 tmp2 tmp3 tmp4
