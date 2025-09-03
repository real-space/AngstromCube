# This file is part of AngstromCube under MIT License

### This makefile for ./a43 should only be used if CMake fails!

HOSTNAME := $(shell hostname)

# NO_UNIT_TESTS=1

# For reproducability compile the git key (hash) into the code
GITKEY=$(shell git log | head -1 | sed -e 's/commit //g')

EXECUTABLE = a43
all: ${EXECUTABLE}

OBJDIR=obj_a43

OBJATOM = $(addprefix $(OBJDIR)/, \
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
	unit_system.o \
	single_atom.o )

OBJ = $(addprefix $(OBJDIR)/, \
	action_plan.o \
	atom_communication.o \
	dyadic_plan.o \
	kinetic_plan.o \
	geometry_input.o \
	geometry_analysis.o \
	finite_difference.o \
	self_consistency.o \
	iterative_poisson.o \
	inline_math.o \
	parallel_poisson.o \
	bessel_transform.o \
	dense_solver.o \
	multi_grid.o \
	mpi_parallel.o \
	potential_generator.o \
	parallel_potential.o \
	poisson_solver.o \
	plane_wave.o \
	sho_basis.o \
	sho_overlap.o \
	sho_hamiltonian.o \
	simple_timer.o \
	energy_mesh.o \
	energy_contour.o \
	sho_projection.o \
	sho_potential.o \
	structure_solver.o \
	green_parallel.o \
	green_action.o \
	green_input.o \
	green_kinetic.o \
	green_dyadic.o \
	green_experiments.o \
	green_sparse.o \
	green_solver.o \
	green_function.o \
	green_memory.o \
	green_potential.o \
	green_tests.o \
	load_balancer.o \
	a43_main.o)

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
INC += -I../external/rapidxml
INC += -I../external/rapidjson

ERRFLAGS =
FEAFLAGS =
FEAFLAGS += -D_GIT_KEY=$(GITKEY)
# FEAFLAGS += -DNO_UNIT_TESTS
# FEAFLAGS += -D_Output_Units_Fixed

# NDEBUG will disable all the assert() checks from <cassert> 
# FEAFLAGS += -DNDEBUG

FEAFLAGS += -DHAS_NO_MPI
FEAFLAGS += -DHAS_NO_OMP
FEAFLAGS += -DHAS_NO_CUDA
FEAFLAGS += -DHAS_BITMAP_EXPORT
FEAFLAGS += -DHAS_SINGLE_ATOM
FEAFLAGS += -DENERGY_CONTOUR_SVG_EXPORT
FEAFLAGS += -DGREEN_FUNCTION_SVG_EXPORT
# FEAFLAGS += -DUSE_SIMPLE_RANGE_TRUNCATION
# FEAFLAGS += -DUSE_RECIPROCAL_RADIAL_GRID
# FEAFLAGS += -DLARGE_ANGULAR_GRIDS
FEAFLAGS += -DGENERAL_CELL
FEAFLAGS += -DHAS_TFQMRGPU
INC += -I/Users/pbaum/tfQMRgpu/include

FEAFLAGS += -DHAS_RAPIDXML
FEAFLAGS += -DHAS_RAPIDJSON
# FEAFLAGS += -DHAS_NO_LAPACK
# FEAFLAGS += -DHAS_NO_BLAS

# FEAFLAGS += -DDEBUGGPU
# ERRFLAGS = -Werror -Wall -Wfatal-errors
# FEAFLAGS += -fopenmp -fno-omit-frame-pointer

ifdef NO_UNIT_TESTS
	## for production: eliminate the unit tests
	FEAFLAGS += -DNO_UNIT_TESTS
endif

## include also code which is not meant for release,
## use github.com/BR903/cppp -UDEVEL to do partial preprocessing
FEAFLAGS += -DDEVEL

OPTFLAGS=

# OPTFLAGS += -Ofast
# OPTFLAGS += -O2
OPTFLAGS += -O0

#### code coverage: compiling generates one *.gcno file per object, exectution generates one *.gcda file per object
####                running gcov generates a *.gcov file per source and per relevant header file
####                see https://danielsieger.com/blog/2022/03/06/code-coverage-for-cpp.html
# OPTFLAGS +=  -O0 -coverage

OPTFLAGS += -g -pedantic -Wall -Wextra -Wno-unused-parameter
# gcc-options
OPTFLAGS += -Wno-format
# OPTFLAGS += -Wno-format-security
OPTFLAGS += -Wno-sign-compare

# for Clang: AddressSanitizer
OPTFLAGS += -fsanitize=address

# for Clang: UndefinedBehaviorSanitizer
OPTFLAGS += -fsanitize=undefined

# LDFLAGS = -mkl
LDFLAGS = ## -lmkl_intel_thread


EXTFLAGS ?=

# ifneq ( ,$(findstring MacBook-Pro,$(HOSTNAME)))
EXTFLAGS += -DHAS_NO_MKL
LDFLAGS += -lm
LDFLAGS += -lblas
# OPTFLAGS += -fopenmp ### works only with gcc
LDFLAGS += -llapack
EXTFLAGS += -DHAS_FFTW
# INC += -I/opt/homebrew/Cellar/fftw/3.3.10_1/include  ### needed by gcc
LDFLAGS += -lfftw3 -L/opt/homebrew/Cellar/fftw/3.3.10_1/lib
#endif

	# on JUSUF:   module --force purge && module load Stages/2022 GCC/11.2.0 OpenMPI/4.1.2 imkl/2021.4.0
	# LDFLAGS += -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

FLAGS = ${ERRFLAGS} ${FEAFLAGS} ${OPTFLAGS} ${EXTFLAGS}

CXXFLAGS = -std=c++11 ${FLAGS}

# the collection of all dependencies
DEP=$(OBJ:%.o=%.d) $(OBJATOM:%.o=%.d)

# target to build executable a43: requirements *.o, envokes the linker
${EXECUTABLE}: ${OBJATOM} ${OBJ}
		${CXX} ${CXXFLAGS} $^ ${LDFLAGS} -o $@

# target to build the static library libliveatom.a 
libliveatom.a: ${OBJATOM} $(OBJDIR)/library_kind_static.o
		ar rcsv $@ $^

# target to build the shared object library libliveatom.so
libliveatom.so: ${OBJATOM} $(OBJDIR)/library_kind_dynamic.o
		${CXX} ${CXXFLAGS} -fPIC -shared $^ ${LDFLAGS} -o $@

# rule how to compile a *.cxx-source file into an *.o-object file
$(OBJDIR)/%.o:%.cxx $(OBJDIR)/%.d
		${CXX} ${CXXFLAGS} ${INC} -c $< -o $@

# rule to generate *.d-dependency files, uses the compiler with all preprocessor flags
$(OBJDIR)/%.d:%.cxx
		${CXX} ${CXXFLAGS} ${INC} -MM $^ -MT $(patsubst %.d,%.o,$@) -MF $@ 

# include dependencies here to know which parts need to be recompiled on change
-include ${DEP}
