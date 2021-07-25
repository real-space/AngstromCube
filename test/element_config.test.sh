#!/usr/bin/env bash

(cd ../src/ && make -j)
exe=../src/a43

## make sure A43's single_atom.cxx has been compiled with -D HAS_AUTO_CONF

## good semicore examples are Na (Z=11, 2s), Zn (Z=30, 3p) and ue (Z=119, 6d)

## experiment: check the semicore feature with sodium Z=11 Na
$exe  +verbosity=7 \
    -t single_atom \
      +single_atom.test.Z=$1 \
      +single_atom.config=element_config \
      +single_atom.core.state.localization=-1 \
      +element_config.core.semicore=-4.5 \
      +element_config.semicore.valence=-1.5 \
      +element_config.rcut=2.26 \
      +element_config.sigma=0.7 \
      +element_config.numax=2 \
      +control.show=1 \
    >  element_config.semicore.$1.out

## We can use the energy criterion to separate core and valence states
## with an energy window of 0 for semicore states by using this:
#     +element_config.core.valence=-2.2 \

$exe  +verbosity=7 \
    -t single_atom \
      +single_atom.test.Z=11 \
      +single_atom.config=sigma_config \
      +element_Na="3s* 1 2p 6 2sSemicore 2 3d | 2.26 sigma 0.7" \
      +control.show=1 \
    >  element_config.semicore.Na.out

## if 11 == $1    diff element_config.semicore.11.out element_config.semicore.Na.out


# exit
# now test a 1s core hole in gold

$exe  +verbosity=7 \
    -t single_atom \
      +single_atom.test.Z=79 \
      +single_atom.config=element_config \
      +single_atom.core.state.localization=-1 \
      +single_atom.test.ion=.5 \
      +element_config.core.valence=-1.5 \
      +element_config.rcut=2.5 \
      +element_config.sigma=.667 \
      +element_config.numax=2 \
      +element_Au.hole.enn=1 \
      +element_Au.hole.ell=0 \
      +element_Au.hole.charge=.5 \
      +control.show=1 \
    >  element_config.corehole.79.out

# default:
#     +element_Au="6s* 1 0 6p* 2e-99 5d* 10 | 2.5 sigma .667"
$exe  +verbosity=7 \
    -t single_atom \
      +single_atom.test.Z=79 \
      +single_atom.config=sigma_config \
      +element_Au="1sCore 1 .5 6s* 1 .5 6p 5d 10 | 2.5 sigma .667" \
      +control.show=1 \
    >  element_config.corehole.Au.out
