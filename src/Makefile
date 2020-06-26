SYSTEM ?= zam044

## Check whether CUDA compiler is available
ENABLEGPU := $(shell command -v nvcc 2> /dev/null)
#ENABLEGPU := $() # Use this line to test non-CUDA version on a system with CUDA installed

# For reproducability compile the git key (hash) into the code
GITKEY=$(shell git log | head -1 | sed -e 's/commit //g')

THISMACHINE := $(shell hostname)

EXECUTABLE = a43
all: ${EXECUTABLE}

OBJ = \
	radial_grid.o \
	radial_integrator.o \
	radial_eigensolver.o \
	radial_potential.o \
	exchange_correlation.o \
	conjugate_gradients.o \
	geometry_analysis.o \
	bessel_transform.o \
	iterative_poisson.o \
	real_space.o \
	davidson_solver.o \
	grid_operators.o \
	finite_difference.o \
	recorded_warnings.o \
	potential_generator.o \
	density_generator.o \
	chemical_symbol.o \
	sigma_config.o \
	fourier_poisson.o \
	atom_core.o \
	single_atom.o \
	sho_overlap.o \
	sho_potential.o \
	sho_projection.o \
	sho_radial.o \
	sho_tools.o \
	sho_unitary.o \
	angular_grid.o \
	control.o \
	main.o

ifdef ENABLEGPU
	OBJ += gpu_drivers.o
endif


.PHONY: clean allclean
clean:
	-rm -f *.o *~ *.d

allclean: clean
	-rm -f ${EXECUTABLE}

CU     = nvcc
### -G
CC     = gcc
# CXX    = icpc
CXX    = g++
##MPI
# CC  = mpicc
# CXX = mpic++

SRC ?= ./src
INC = -I/usr/local/cuda/include


ERRFLAGS =
FEAFLAGS =
# FEAFLAGS += -D_Output_Units_Fixed
FEAFLAGS += -D_GIT_KEY=$(GITKEY) 
# ERRFLAGS = -Werror -Wall -Wfatal-errors
FEAFLAGS += ###-fopenmp -fno-omit-frame-pointer
# FEAFLAGS += -D NO_UNIT_TESTS ## for production

## include also code which is not meant for release,
## use github.com/BR903/cppp -UDEVEL to do partial preprocessing
FEAFLAGS += -D DEVEL 

# normal mode
# OPTFLAGS = -Ofast 
OPTFLAGS = -O0
### -mtune=power8 -mcpu=power8 -mpower8-fusion -mpower8-vector -mvsx -maltivec -mdirect-move -fstrict-aliasing -fomit-frame-pointer -fno-schedule-insns
# OPTFLAGS = -O2
OPTFLAGS += -g -pedantic -Wall -Warray-bounds -Wno-sign-compare -Wno-nested-anon-types -Wno-format-security -Wno-format

# LDFLAGS = -mkl
LDFLAGS = ## -lmkl_intel_thread

CUOPTFLAGS =
#### production mode
CUOPTFLAGS += -O3 --use_fast_math
#### debug mode
# OPTFLAGS  = -O0 -g
# CUOPTFLAGS += -g -O0 -G
# CUOPTFLAGS += -DDEBUGGPU

## verbose assembler
CUOPTFLAGS +=  -Xptxas -v

# for profiling
CUOPTFLAGS += -lineinfo
ifdef ENABLEGPU
	FEAFLAGS += -DUSE_NVTX
	LDFLAGS  += -rdynamic -L/usr/local/cuda/lib64 -lcudart -lnvToolsExt
endif

EXTFLAGS ?=

ifeq ($(SYSTEM),zam044)
	EXTFLAGS += -D HAS_no_MKL
	LDFLAGS += -lm -lblas -llapack ### works on zam044 (MacBook Pro)
else
	## eq ($(THISMACHINE), "zam054")
	### works with g++ on zam054 (JSC workstation)
	LDFLAGS += -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread 
endif

FLAGS = ${ERRFLAGS} ${FEAFLAGS} ${OPTFLAGS} ${EXTFLAGS}

GENCODE_SM60= -gencode arch=compute_60,code=sm_60 -gencode arch=compute_60,code=compute_60
GENCODE_SM70= -gencode arch=compute_70,code=sm_70 -gencode arch=compute_70,code=compute_70
GENCODE_FLAGS= ${GENCODE_SM60} ${GENCODE_SM70}

CCFLAGS  = -std=gnu99 ${FLAGS}
# K80
#CUFLAGS  = -std=c++11 ${GENCODE_FLAGS} ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
# P100
CUFLAGS  = -std=c++11 ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
CUFLAGS  = -std=c++11 ${GENCODE_FLAGS} ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
CXXFLAGS = -std=c++11 ${FLAGS}


DEP=$(OBJ:%.o=%.d)

${EXECUTABLE}: ${OBJ}
	${CXX} ${CXXFLAGS} $^ ${LDFLAGS}  -o $@

%.o:${SRC}/%.c %.d
	${CC} ${CCFLAGS} ${INC} -c $<
%.o:${SRC}/%.cu %.d
	${CU} ${CUFLAGS} ${INC} -c $<
%.o:${SRC}/%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<

%.o:%.c %.d
	${CC} ${CCFLAGS} ${INC} -c $<
%.o:%.cu %.d
	${CU} ${CUFLAGS} ${INC} -c $<
%.o:%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<

%.d:${SRC}/%.c
	${CC} ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:${SRC}/%.cu
	${CC} -x c++ ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:${SRC}/%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@

%.d:%.c
	${CC} ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cu
	${CC} -x c++ ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@

-include ${DEP}
