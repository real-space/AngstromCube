SYSTEM ?= zam044

## Check whether CUDA compiler is available
ENABLE_GPU := $(shell command -v nvcc 2> /dev/null)
#ENABLE_GPU := $() # Use this line to test non-CUDA version on a system with CUDA installed

# NO_UNIT_TESTS=1

# For reproducability compile the git key (hash) into the code
GITKEY=$(shell git log | head -1 | sed -e 's/commit //g')

THISMACHINE := $(shell hostname)

EXECUTABLE = a43
all: ${EXECUTABLE}

OBJATOM = \
	control.o \
	recorded_warnings.o \
	radial_grid.o \
	radial_integrator.o \
	radial_potential.o \
	radial_eigensolver.o \
	exchange_correlation.o \
	chemical_symbol.o \
	display_units.o \
	sigma_config.o \
	atom_core.o \
	sho_unitary.o \
	angular_grid.o \
	single_atom.o

OBJ = \
	geometry_analysis.o \
	self_consistency.o \
	iterative_poisson.o \
	bessel_transform.o \
	potential_generator.o \
	poisson_solver.o \
	plane_wave.o \
	sho_basis.o \
	sho_overlap.o \
	sho_hamiltonian.o \
	sho_projection.o \
	structure_solver.o \
	main.o


ifdef ENABLE_GPU
	OBJ += \
	gpu_drivers.o
endif


.PHONY: clean allclean
clean:
	-rm -f *.o *~ *.d *.a

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

INC ?=
INC += -I../include
INC += -I../external
# INC += -I/usr/local/cuda/include


ERRFLAGS =
FEAFLAGS =
# FEAFLAGS += -D NO_UNIT_TESTS
# FEAFLAGS += -D_Output_Units_Fixed
FEAFLAGS += -D_GIT_KEY=$(GITKEY)
FEAFLAGS += -D HAS_NO_MPI
FEAFLAGS += -D HAS_NO_CUDA
# FEAFLAGS += -D USE_RECIPROCAL_RADIAL_GRID
# FEAFLAGS += -D LARGE_ANGULAR_GRIDS
FEAFLAGS += -D HAS_RAPIDXML
FEAFLAGS += -D HAS_RAPIDJSON
# ERRFLAGS = -Werror -Wall -Wfatal-errors
# FEAFLAGS += -fopenmp -fno-omit-frame-pointer

ifdef NO_UNIT_TESTS
	## for production: eliminate the unit tests
	FEAFLAGS += -D NO_UNIT_TESTS
endif

## include also code which is not meant for release,
## use github.com/BR903/cppp -UDEVEL to do partial preprocessing
FEAFLAGS += -D DEVEL

# normal mode
# OPTFLAGS = -Ofast
# OPTFLAGS = -O2
OPTFLAGS = -O0
OPTFLAGS += -g -pedantic -Wall
### -mtune=power8 -mcpu=power8 -mpower8-fusion -mpower8-vector -mvsx -maltivec -mdirect-move -fstrict-aliasing -fomit-frame-pointer -fno-schedule-insns
# -Wno-format-security
# -Wno-format

# LDFLAGS = -mkl
LDFLAGS = ## -lmkl_intel_thread

CUOPTFLAGS =
#### production mode
CUOPTFLAGS += -O3 --use_fast_math
#### debug mode
# CUOPTFLAGS += -g -O0 -G
# CUOPTFLAGS += -DDEBUGGPU

## verbose assembler
CUOPTFLAGS +=  -Xptxas -v

# for profiling
CUOPTFLAGS += -lineinfo
ifdef ENABLE_GPU
	FEAFLAGS += -DUSE_NVTX
	LDFLAGS  += -rdynamic -L/usr/local/cuda/lib64 -lcudart -lnvToolsExt
endif

EXTFLAGS ?=

ifeq ($(SYSTEM),zam044)
	EXTFLAGS += -D HAS_NO_MKL
	EXTFLAGS += -D HAS_FFTW
	LDFLAGS += -lm -lblas -llapack ### works on zam044 (MacBook Pro)
	LDFLAGS += -lfftw3
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


DEP=$(OBJ:%.o=%.d) $(OBJATOM:%.o=%.d)

${EXECUTABLE}: ${OBJATOM} ${OBJ}
	${CXX} ${CXXFLAGS} $^ ${LDFLAGS}  -o $@

libliveatom.a: ${OBJATOM}
	ar rcsv $@ $^

libliveatom.so: ${OBJATOM}
	${CXX} ${CXXFLAGS} -fPIC -shared $^ ${LDFLAGS} -o $@

%.o:%.c %.d
	${CC} ${CCFLAGS} ${INC} -c $<
%.o:%.cu %.d
	${CU} ${CUFLAGS} ${INC} -c $<
%.o:%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<

%.d:%.c
	${CC} ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cu
	${CC} -x c++ ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@

-include ${DEP}
