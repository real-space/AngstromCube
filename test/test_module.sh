#!/usr/bin/env bash

module=$1

sed -e "s/MODULE/$module/g" test_MODULE.cxx > test_$module.cxx
g++ -std=c++11 \
    -I../include \
    -O0 \
    -g \
    -pedantic \
    -Wall \
    -Wno-format-security \
    -Wno-format \
    -D STANDALONE_TEST \
    test_$module.cxx \
    ../src/recorded_warnings.o \
    ../src/control.o \
    -o test_$module.exe \
 && ./test_$module.exe $@

rm  -f test_$module.cxx test_$module.exe 
rm -rf test_$module.exe.dSYM
