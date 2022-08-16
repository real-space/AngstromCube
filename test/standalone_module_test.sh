#!/usr/bin/env bash

### This script tests if a header file <module>.hxx has all the
### include statements it needs to be compiled standalone.
### Linking will require a dependency tree of objects, so we skip that.
###
### This script is called by standalone_module_tests.sh

    module=$1
    echo $module

    ## generate a short main function including nothing but the header
    echo "#include \"$module.hxx\"    // ::all_tests"         > test_me.cxx
    echo "int main() { return int($module::all_tests(6)); }" >> test_me.cxx

    ## compile to check for missing include files
    g++ -std=c++11 \
        -I../include/ \
        -I../external/ \
        -g -pedantic -Wall -O0 \
           -Wno-sign-compare \
           -Wno-format \
              -D HAS_NO_MKL \
              -D DEVEL \
              -D HAS_NO_MPI \
              -D HAS_RAPIDXML \
              -D HAS_NO_CUDA \
              -c test_me.cxx

#               -D NO_UNIT_TESTS \

    ## cleanup
    rm -f test_me.cxx test_me.o
