# This file is part of AngstromCube under MIT License


### JUSUF:   module --force purge && module load Stages/2020 GCC/10.3.0 OpenMPI/4.1.1 CUDA/11.3

## Check whether CUDA compiler is available
ENABLE_GPU := $(shell command -v nvcc 2> /dev/null)
#ENABLE_GPU := $() # Use this line to test non-CUDA version on a system with CUDA installed

# For reproducability compile the git key (hash) into the code
GITKEY=$(shell git log | head -1 | sed -e 's/commit //g')

HOSTNAME := $(shell hostname)

EXECUTABLE = green
all: ${EXECUTABLE}

OBJECTS = \
	control.o \
	recorded_warnings.o \
	display_units.o \
	load_balancer.o \
	mpi_parallel.o \
	green_experiments.o \
	green_function.o \
	green_parallel.o \
	green_kinetic.o \
	green_input.o \
	green_tests.o \
	green.o


.PHONY: clean allclean
clean:
	-rm -f *.o *~ *.d

allclean: clean
	-rm -f ${EXECUTABLE}

CU    ?= nvcc
CC    ?= gcc
FEAFLAGS ?=
INC ?=
LDFLAGS ?=
ERRFLAGS ?=
OPTFLAGS ?=
CUOPTFLAGS ?=
CXXOPTFLAGS ?=

ifeq ($(HOSTNAME),MacBook-Pro.local)
	### Apple clang version 13.1.6 (clang-1316.0.21.2.5)
	CXX    = g++ -Wno-unused-command-line-argument
	FEAFLAGS += -DHAS_NO_MPI
	FEAFLAGS += -DHAS_LAPACK
	EXTFLAGS += -DHAS_NO_MKL
	LDFLAGS  += -lm -lblas -llapack
else
	### JUSUF: module --force purge && module load Stages/2020 GCC/10.3.0 OpenMPI/4.1.1 CUDA/11.3
	CXX    = mpic++ -Wno-unused-command-line-argument
	FEAFLAGS += -DHAS_NO_MKL
	LDFLAGS  += -lmpi
endif


INC += -I../include
INC += -I../external


# FEAFLAGS += -D_Output_Units_Fixed
FEAFLAGS += -D_GIT_KEY=$(GITKEY)
# FEAFLAGS += -DHAS_NO_MPI

FEAFLAGS += -DHAS_RAPIDXML
FEAFLAGS += -DHAS_RAPIDJSON
FEAFLAGS += -DUSE_SIMPLE_RANGE_TRUNCATION
# FEAFLAGS += -DNO_UNIT_TESTS

## use the tfQMRgpu library to invert the shifted Hamiltonian
FEAFLAGS += -DHAS_TFQMRGPU
INC += -I../../tfQMRgpu

ERRFLAGS =
# ERRFLAGS += -Werror
ERRFLAGS += -Wall -Wfatal-errors
# FEAFLAGS += -fopenmp -fno-omit-frame-pointer

## for production: eliminate the unit tests
# FEAFLAGS += -DNO_UNIT_TESTS

## include also code which is not meant for release,
## use github.com/BR903/cppp -UDEVEL to do partial preprocessing
FEAFLAGS += -DDEVEL


# normal mode
# OPTFLAGS += -Ofast
# OPTFLAGS += -O2
# OPTFLAGS += -O0
OPTFLAGS += -g -Wall
CXXOPTFLAGS += -pedantic
# gcc-options
# OPTFLAGS += -Wno-format-security
OPTFLAGS += -Wno-sign-compare
OPTFLAGS += -Wno-format
OPTFLAGS += -Wno-unknown-pragmas

# LDFLAGS += -mkl
LDFLAGS += ## -lmkl_intel_thread

CUOPTFLAGS =
#### production mode
# CUOPTFLAGS += -O3 --use_fast_math
CUOPTFLAGS += -O2
#### debug mode
# CUOPTFLAGS += -g -O0 -G -lineinfo
# CUOPTFLAGS += -DDEBUGGPU

## verbose GPU assembler
CUOPTFLAGS += -Xptxas -v
CUOPTFLAGS += --expt-relaxed-constexpr

# for profiling
CUOPTFLAGS += -lineinfo
ifdef ENABLE_GPU
	INC += -I/usr/local/cuda/include
	FEAFLAGS += -DUSE_NVTX
	LDFLAGS  += -L/usr/local/cuda/lib64 -lcudart -lnvToolsExt -lcurand #-rdynamic
#	FEAFLAGS += -DHAS_NO_CUDA
else
	FEAFLAGS += -DHAS_NO_CUDA
endif

EXTFLAGS ?=

FLAGS = ${ERRFLAGS} ${FEAFLAGS} ${OPTFLAGS} ${EXTFLAGS}

GENCODE60= -gencode arch=compute_60,code=sm_60 -gencode arch=compute_60,code=compute_60 ## for P100
GENCODE70= -gencode arch=compute_70,code=sm_70 -gencode arch=compute_70,code=compute_70 ## for V100
GENCODE80= -gencode arch=compute_80,code=sm_80 -gencode arch=compute_80,code=compute_80 ## for A100
GENCODE_FLAGS=
# GENCODE_FLAGS += ${GENCODE60}
GENCODE_FLAGS += ${GENCODE70}
# GENCODE_FLAGS += ${GENCODE80}

CCFLAGS  = -std=gnu99 ${FLAGS}
# K80
#CUFLAGS  = -std=c++11 ${GENCODE_FLAGS} ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
# P100
CUFLAGS  = -std=c++11 ${GENCODE_FLAGS} ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
CXXFLAGS = -std=c++11 ${FLAGS} ${CXXOPTFLAGS}

ifdef ENABLE_GPU
else
	CU = $(CXX)
	CUFLAGS = $(CXXFLAGS)
endif


DEP=$(OBJECTS:%.o=%.d)

${EXECUTABLE}: ${OBJECTS}
	${CU} ${CUFLAGS} $^ ${LDFLAGS}  -o $@

%.o:%.c %.d
	${CC} ${CCFLAGS} ${INC} -c $<
%.o:%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<
%.o:%.cu %.d
	${CU} ${CUFLAGS} ${INC} -c $<

%.d:%.c
	${CC} ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cu
	${CU} ${CUFLAGS} ${INC} -MM $^ -MF $@

-include ${DEP}
