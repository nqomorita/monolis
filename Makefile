#> monolis_utils Makefile

##> compiler setting
FC     = mpif90
FFLAGS = -fPIC -O2 -mtune=native -std=legacy -Wno-missing-include-dirs
CC     = mpicc -std=c99
CFLAGS = -fPIC -O2
LINK   = $(FC)

##> directory setting
MOD_DIR = -J ./include
BIN_DIR = ./bin
SRC_DIR = ./src
TST_DIR = ./src_test
OBJ_DIR = ./obj
LIB_DIR = ./lib
WRAP_DIR= ./wrapper
TST_WRAP_DIR = ./wrapper_test
DRV_DIR = ./driver
LIBRARY = libmonolis.a
LIBRARY_SOLVER = libmonolis_solver.a
CPP     = -cpp
CPPFLAG =

#INCLUDE = -I /Users/morita/opt/include -I ./include -I /usr/include -I ./submodule/gedatsu/include -I ./submodule/monolis_utils/include
INCLUDE = -I ./include -I /usr/include -I ./submodule/gedatsu/include -I ./submodule/monolis_utils/include

##> compiler option setting
ifdef FLAGS
	comma:= ,
	empty:=
	space:= $(empty) $(empty)
	DFLAGS = $(subst $(comma), $(space), $(FLAGS))

	ifeq ($(findstring DEBUG, $(DFLAGS)), DEBUG)
		FFLAGS  = -fPIC -O2 -std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow -Wno-missing-include-dirs
		CFLAGS  = -fPIC -O2 -g -ggdb
	endif

	ifeq ($(findstring INTEL, $(DFLAGS)), INTEL)
		FC      = mpiifx -qmkl=cluster
		FFLAGS  = -fPIC -O2 -align array64byte -nofor-main
		CC      = mpiicx
		CFLAGS  = -fPIC -O2 -no-multibyte-chars
		MOD_DIR = -module ./include
		LINK    = $(FC)
	endif

	ifeq ($(findstring A64FX, $(DFLAGS)), A64FX)
		FC      = mpifrtpx
		FFLAGS  = -Nalloc_assign -Kfast -SCALAPACK -SSL2
		CC      = mpifccpx -Nclang 
		CFLAGS  = -Kfast
		MOD_DIR = -M ./include
		LINK    = mpiFCCpx --linkfortran -SSL2
		#INCLUDE = -I ./include -I ./submodule/gedatsu/include -I ./submodule/monolis_utils/include
	endif
endif

USE_LIB_CORE = -L./lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis
#USE_LIB_OPT  = -L/Users/morita/opt/lib -lscalapack -lopenblas -lc++
USE_LIB_OPT  = -L./lib -lscalapack -llapack -lblas

##> liblary option setting
ifdef FLAGS
	comma:= ,
	empty:=
	space:= $(empty) $(empty)
	DFLAGS = $(subst $(comma), $(space), $(FLAGS))

	ifeq ($(findstring ML, $(DFLAGS)), ML)
		CPPFLAG += -DWITH_ML
		#INCLUDE_ML = 
		USE_LIB_ML = -L./lib -lml -lzoltan -lmetis
	endif

	ifeq ($(findstring MUMPS, $(DFLAGS)), MUMPS)
		CPPFLAG += -DWITH_MUMPS
		#INCLUDE_MUMPS = 
		USE_LIB_MUMPS = -L./lib -ldmumps -lmumps_common -lpord
	endif

	ifeq ($(findstring BLOPEX, $(DFLAGS)), BLOPEX)
		#INCLUDE_BLOPEX = 
		CPPFLAG += -DWITH_BLOPEX
		USE_LIB_BLOPEX = -L./lib 
	endif
endif

USE_LIB = $(USE_LIB_CORE) $(USE_LIB_ML) $(USE_LIB_MUMPS) $(USE_LIB_BLOPEX) $(USE_LIB_OPT)

##> other commands
MAKE = make
CD   = cd
CP   = cp
RM   = rm -rf
AR   = - ar ruv

##> **********
##> target (1)
LIB_TARGET = $(LIB_DIR)/$(LIBRARY_SOLVER)

##> source file define
SRC_DEFINE = \
def_solver_prm.f90 \
def_mat.f90 \
def_struc.f90 \
def_solver_prm_util.f90

SRC_MAT = \
spmat_copy.f90 \
spmat_handler_util.f90 \
spmat_nzpattern_util.f90 \
spmat_reordering.f90 \
spmat_nzpattern.f90 \
spmat_handler.f90 \
spmat_convert_sym.f90

SRC_LINALG = \
inner_product.f90 \
vec_util.f90 \
matvec.f90 \
mat_converge.f90 

SRC_WRAP = \
wrapper_lapack.f90 \
wrapper_scalapack.f90 \
monolis_wrapper_ml.c \
monolis_wrapper_blopex.c 

SRC_FACT = \
fillin.f90 \
analysis.f90 \
factorize.f90 \
LU/LU_nn.f90

SRC_OPT = \
nnls.f90 

SRC_PREC = \
diag/diag_33.f90 \
diag/diag_nn.f90 \
diag/diag_V.f90 \
sor/sor_33.f90 \
sor/sor_nn.f90 \
sor/sor_V.f90 \
diag.f90 \
sor.f90 \
LU.f90 \
MUMPS.f90 \
ML.f90 \
precond.f90

SRC_ITER = \
CG.f90 \
BiCGSTAB.f90 \
BiCGSTAB_N128.f90 \
BiCGSTAB_noprec.f90 \
GropCG.f90 \
PipeCG.f90 \
PipeCR.f90 \
PipeBiCGSTAB.f90 \
PipeBiCGSTAB_noprec.f90 \
SOR.f90 \
BiCGSAFE.f90 \
IDRS.f90 \
COCG.f90

SRC_SOLV = \
solver.f90

SRC_EIGEN = \
Lanczos_util.f90 \
Lanczos.f90 \
LOBPCG_util.f90 \
LOBPCG.f90 \
eigen_solver.f90

SRC_ALL = \
$(addprefix define/, $(SRC_DEFINE)) \
$(addprefix matrix/, $(SRC_MAT)) \
$(addprefix linalg/, $(SRC_LINALG)) \
$(addprefix wrapper/, $(SRC_WRAP)) \
$(addprefix fact/, $(SRC_FACT)) \
$(addprefix prec/, $(SRC_PREC)) \
$(addprefix iterative/, $(SRC_ITER)) \
$(addprefix solver/, $(SRC_SOLV)) \
$(addprefix eigen/, $(SRC_EIGEN)) \
$(addprefix optimize/, $(SRC_OPT))

##> C wrapper section
SRC_DEFINE_C = \
monolis_def_solver_prm_c.c \
monolis_def_mat_c.c \
monolis_def_struc_c.c \
monolis_def_solver_prm_util_c.c

SRC_LINALG_CF = \
matvec_wrap.f90 \
inner_product_wrap.f90

SRC_LINALG_C = \
monolis_matvec_c.c \
monolis_vec_util_c.c \
monolis_inner_product_c.c

SRC_MAT_CF = \
spmat_nzpattern_util_wrap.f90 \
spmat_handler_util_wrap.f90

SRC_MAT_C = \
monolis_spmat_copy_c.c \
monolis_spmat_nzpattern_util_c.c \
monolis_spmat_nzpattern_c.c \
monolis_spmat_handler_c.c \
monolis_spmat_handler_util_c.c

SRC_WRAP_CF = \
scalapack_wrapper.f90 

SRC_WRAP_C = \
monolis_wrapper_scalapack_c.c

SRC_OPT_CF = \
nnls_wrapper.f90 

SRC_OPT_C = \
monolis_nnls_c.c 

SRC_SOLV_CF = \
solver_wrapper.f90

SRC_SOLV_C = \
monolis_solver_c.c 

SRC_EIGEN_SOLV_CF = \
eigen_solver_wrapper.f90

SRC_EIGEN_SOLV_C = \
monolis_eigen_solver_c.c

SRC_ALL_CC = \
$(addprefix define/, $(SRC_DEFINE_C)) \
$(addprefix linalg/, $(SRC_LINALG_C)) \
$(addprefix matrix/, $(SRC_MAT_C)) \
$(addprefix wrapper/, $(SRC_WRAP_C)) \
$(addprefix optimize/, $(SRC_OPT_C)) \
$(addprefix solver/, $(SRC_SOLV_C)) \
$(addprefix eigen/, $(SRC_EIGEN_SOLV_C)) \

SRC_ALL_CF = \
$(addprefix define/, $(SRC_DEFINE_CF)) \
$(addprefix linalg/, $(SRC_LINALG_CF)) \
$(addprefix matrix/, $(SRC_MAT_CF)) \
$(addprefix wrapper/, $(SRC_WRAP_CF)) \
$(addprefix optimize/, $(SRC_OPT_CF)) \
$(addprefix solver/, $(SRC_SOLV_CF)) \
$(addprefix eigen/, $(SRC_EIGEN_SOLV_CF))

SRC_ALL_C = \
$(SRC_ALL_CF) \
$(SRC_ALL_CC)

##> lib objs
LIB_SOURCES = \
$(addprefix $(SRC_DIR)/,  $(SRC_ALL)) \
$(addprefix $(WRAP_DIR)/, $(SRC_ALL_C)) \
./src/monolis_solver.f90 \
./src/monolis.f90
LIB_OBJS1   = $(subst $(SRC_DIR), $(OBJ_DIR), $(LIB_SOURCES:.f90=.o))
LIB_OBJS    = $(subst $(WRAP_DIR), $(OBJ_DIR), $(LIB_OBJS1:.c=.o))

##> **********
##> test target (2)
TEST_TARGET = $(TST_DIR)/monolis_test

##> lib objs
TST_SRC_ALL = $(SRC_ALL) monolis.f90
TST_SOURCES = $(addprefix $(TST_DIR)/, $(TST_SRC_ALL))
TST_OBJSt   = $(subst $(TST_DIR), $(OBJ_DIR), $(TST_SOURCES:.f90=_test.o))
TST_OBJS    = $(TST_OBJSt:.c=_test.o)

##> **********
##> test target (3)
TEST_C_TARGET = $(TST_WRAP_DIR)/monolis_c_test

##> lib objs
SRC_DEFINE_C_TEST = \
monolis_def_solver_prm_c_test.c \
monolis_def_solver_prm_util_c_test.c

SRC_LINALG_C_TEST = \
monolis_inner_product_c_test.c \
monolis_vec_util_c_test.c \
monolis_matvec_c_test.c

SRC_MAT_C_TEST = \
monolis_spmat_copy_c_test.c\
monolis_spmat_nzpattern_util_c_test.c\
monolis_spmat_nzpattern_c_test.c \
monolis_spmat_handler_c_test.c \
monolis_spmat_handler_util_c_test.c

SRC_WRAP_C_TEST = \
monolis_wrapper_scalapack_c_test.c

SRC_NNLS_C_TEST = \
monolis_nnls_c_test.c

SRC_SOLVER_C_TEST = \
monolis_solver_c_test.c

SRC_EIGEN_SOLVER_C_TEST = \
monolis_eigen_solver_c_test.c

SRC_ALL_C_TEST = \
$(addprefix define/, $(SRC_DEFINE_C_TEST)) \
$(addprefix linalg/, $(SRC_LINALG_C_TEST)) \
$(addprefix matrix/, $(SRC_MAT_C_TEST)) \
$(addprefix wrapper/, $(SRC_WRAP_C_TEST)) \
$(addprefix optimize/, $(SRC_NNLS_C_TEST)) \
$(addprefix solver/, $(SRC_SOLVER_C_TEST)) \
$(addprefix eigen/, $(SRC_EIGEN_SOLVER_C_TEST))

TST_SRC_C_ALL = $(SRC_ALL_C_TEST) monolis_c_test.c
TST_C_SOURCES = $(addprefix $(TST_WRAP_DIR)/, $(TST_SRC_C_ALL))
TST_C_OBJS    = $(subst $(TST_WRAP_DIR), $(OBJ_DIR), $(TST_C_SOURCES:.c=.o))

##> header section
C_HEADER = \
$(SRC_ALL_CC:.c=.h) \
monolis_solver.h \
monolis.h 

C_HEADER_FILES = \
$(SRC_DEFINE_C:.c=.h) \
$(SRC_LINALG_C:.c=.h) \
$(SRC_MAT_C:.c=.h) \
$(SRC_WRAP_C:.c=.h) \
$(SRC_OPT_C:.c=.h) \
$(SRC_SOLV_C:.c=.h) \
$(SRC_EIGEN_SOLV_C:.c=.h) \
monolis_solver.h \
monolis.h

##> target
all: \
	cp_header \
	cp_header_lib \
	cp_bin_lib \
	$(LIB_TARGET) \
	$(TEST_TARGET) \
	$(TEST_C_TARGET)

lib: \
	cp_header \
	cp_header_lib \
	$(LIB_TARGET)

$(LIB_TARGET): $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS) $(ARC_LIB)

$(TEST_TARGET): $(TST_OBJS)
	$(LINK) $(FFLAGS) $(CPP) $(CPPFLAG) $(INCLUDE) -o $@ $(TST_OBJS) $(USE_LIB)

$(TEST_C_TARGET): $(TST_C_OBJS)
	$(LINK) $(FFLAGS) $(INCLUDE) -o $@ $(TST_C_OBJS) $(USE_LIB)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(CPPFLAG) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(TST_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(CPPFLAG) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(WRAP_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(CPPFLAG) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(CPPFLAG) $(INCLUDE) -o $@ -c $<

$(OBJ_DIR)/%.o: $(WRAP_DIR)/%.c
	$(CC) $(CFLAGS) $(CPPFLAG) $(INCLUDE) -o $@ -c $<

$(OBJ_DIR)/%.o: $(TST_WRAP_DIR)/%.c
	$(CC) $(CFLAGS) $(CPPFLAG) $(INCLUDE) -o $@ -c $<

cp_header:
	$(CP) $(addprefix $(WRAP_DIR)/, $(C_HEADER)) ./include/
	$(CP) ./src/wrapper/monolis_wrapper_ml.h ./include/

cp_header_lib:
	$(CP) ./submodule/monolis_utils/include/* ./include/
	$(CP) ./submodule/monolis_utils/lib/libmonolis_utils.a ./lib/
	$(CP) ./submodule/gedatsu/include/* ./include/
	$(CP) ./submodule/gedatsu/lib/libgedatsu.a ./lib/
	$(CP) ./submodule/ggtools/include/* ./include/
	$(CP) ./submodule/ggtools/lib/libggtools.a ./lib/

cp_bin_lib:
	$(CP) ./submodule/monolis_utils/bin/* ./bin/
	$(CP) ./submodule/gedatsu/bin/* ./bin/
#	$(CP) ./submodule/ggtools/bin/* ./bin/

clean:
	$(RM) $(LIB_OBJS) \
	$(RM) $(TST_OBJS) \
	$(RM) $(TST_C_OBJS) \
	$(RM) $(LIB_TARGET) \
	$(RM) $(TEST_TARGET) \
	$(RM) $(TEST_C_TARGET) \
	$(RM) ./include/*.mod \
	$(RM) $(addprefix ./include/, $(C_HEADER_FILES)) \
	$(RM) ./bin/*

.PHONY: clean
