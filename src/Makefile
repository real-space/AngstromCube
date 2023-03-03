HOSTNAME := $(shell hostname)

# NO_UNIT_TESTS=1

# For reproducability compile the git key (hash) into the code
GITKEY=$(shell git log | head -1 | sed -e 's/commit //g')

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
	scattering_test.o \
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
	dense_solver.o \
	potential_generator.o \
	poisson_solver.o \
	plane_wave.o \
	sho_basis.o \
	sho_overlap.o \
	sho_hamiltonian.o \
	sho_projection.o \
	sho_potential.o \
	structure_solver.o \
	green_parallel.o \
	green_function.o \
	green_experiments.o \
	green_input.o \
	load_balancer.o \
	green_tests.o \
	main.o

.PHONY: clean allclean
clean:
	-rm -f *.o *.so *~ *.d *.a

allclean: clean
	-rm -f ${EXECUTABLE}

CC     = gcc
CXX   ?= g++
##MPI
# CC  = mpicc
# CXX = mpic++

INC ?=
INC += -I../include
INC += -I../external

ERRFLAGS =
FEAFLAGS =
FEAFLAGS += -D _GIT_KEY=$(GITKEY)
FEAFLAGS += -D NO_UNIT_TESTS
# FEAFLAGS += -D _Output_Units_Fixed
FEAFLAGS += -D HAS_NO_MPI
FEAFLAGS += -D HAS_NO_CUDA
FEAFLAGS += -D HAS_BITMAP_EXPORT
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
OPTFLAGS=
# OPTFLAGS += -Ofast
# OPTFLAGS += -O2
OPTFLAGS += -O0
OPTFLAGS += -g -pedantic -Wall
# gcc-options
OPTFLAGS += -Wno-format
# OPTFLAGS += -Wno-format-security
OPTFLAGS += -Wno-sign-compare

# LDFLAGS = -mkl
LDFLAGS = ## -lmkl_intel_thread


EXTFLAGS ?=

ifeq ($(HOSTNAME),MacBook-Pro.local)
	EXTFLAGS += -D HAS_NO_MKL
	LDFLAGS += -lm
	LDFLAGS += -lblas
	# OPTFLAGS += -fopenmp ### works only with gcc
	LDFLAGS += -llapack
	EXTFLAGS += -D HAS_FFTW
  # INC += -I/opt/homebrew/Cellar/fftw/3.3.10_1/include  ### needed by gcc
	LDFLAGS += -lfftw3 -L/opt/homebrew/Cellar/fftw/3.3.10_1/lib
else
	# on JUSUF:   module --force purge && module load Stages/2022 GCC/11.2.0 OpenMPI/4.1.2 imkl/2021.4.0
	LDFLAGS += -lm
	LDFLAGS += -lmkl_intel_lp64
	LDFLAGS += -lmkl_intel_thread
	LDFLAGS += -lmkl_core
	LDFLAGS += -liomp5
	LDFLAGS += -lpthread
endif

FLAGS = ${ERRFLAGS} ${FEAFLAGS} ${OPTFLAGS} ${EXTFLAGS}

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
%.o:%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<

%.d:%.c
	${CC} ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@

-include ${DEP}
