#> monolis_utils Makefile

##> compiler setting
FC     = mpif90
FFLAGS = -fPIC -O2 -mtune=native -march=native -std=legacy -Wno-missing-include-dirs
CC     = mpicc -std=c99
CFLAGS = -fPIC -O2

##> directory setting
MOD_DIR = -J ./include
INCLUDE = -I /usr/include -I ./include -I ./submodule/gedatsu/include -I ./submodule/monolis_utils/include
USE_LIB = -L./lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -lscalapack -llapack -lblas
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
CPP     = -cpp $(FLAG_DEBUG)

##> option setting
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
		FC      = mpiifort -qmkl=cluster
		FFLAGS  = -fPIC -O2 -align array64byte -nofor-main
		CC      = mpiicc
		CFLAGS  = -fPIC -O2 -no-multibyte-chars
		MOD_DIR = -module ./include
		USE_LIB = -L./lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis
	endif
endif

##> other commands
MAKE = make
CD   = cd
CP   = cp
RM   = rm -rf
AR   = - ar ruv

##> **********
##> target (1)
LIB_TARGET = $(LIB_DIR)/$(LIBRARY_SOLVER)
LIBALL_TARGET = $(LIB_DIR)/$(LIBRARY)

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
spmat_nzpattern.f90 \
spmat_handler.f90

#spmat_fillin.f90 \
#spmat_reorder.f90 \
#spmat_scaling.f90 \

SRC_LINALG = \
matvec.f90 \
inner_product.f90 \
mat_converge.f90 \
vec_util.f90

SRC_WRAP = \
wrapper_lapack.f90 \
wrapper_scalapack.f90

#matmat.f90 \

#SRC_FACT = \
#11/fact_LU_11.f90 \
#11/fact_MF_11.f90 \
#33/fact_LU_33.f90 \
#nn/fact_LU_nn.f90 \
#fact_LU.f90 \
#fact_MF.f90

SRC_PREC = \
diag/diag_33.f90 \
diag/diag_nn.f90 \
sor/sor_33.f90 \
sor/sor_nn.f90 \
diag.f90 \
sor.f90 \
precond.f90

#MUMPS.f90 \
#ilu.f90 \
#Jacobi.f90 \
#MF.f90 \

SRC_ITER = \
CG.f90 \
BiCGSTAB.f90 \
BiCGSTAB_noprec.f90 \
GropCG.f90 \
PipeCG.f90 \
PipeCR.f90 \
PipeBiCGSTAB.f90 \
PipeBiCGSTAB_noprec.f90 \
COCG.f90

#CABiCGSTAB_noprec.f90 \
#GMRES.f90 \

SRC_SOLV = \
solver.f90

SRC_EIGEN = \
Lanczos_util.f90 \
Lanczos.f90 \
eigen_solver.f90

##> C wrapper section
SRC_DEFINE_C = \
monolis_def_solver_prm_c.c \
monolis_def_mat_c.c \
monolis_def_struc_c.c \
monolis_def_solver_prm_util_c.c

SRC_LINALG_C = \
matvec_wrap.f90 \
inner_product_wrap.f90 \
monolis_matvec_c.c \
monolis_vec_util_c.c \
monolis_inner_product_c.c

SRC_MAT_C = \
monolis_spmat_copy_c.c \
monolis_spmat_nzpattern_util_c.c \
monolis_spmat_nzpattern_c.c \
monolis_spmat_handler_c.c \
spmat_nzpattern_util_wrap.f90 \
spmat_handler_util_wrap.f90

SRC_WRAP_C = \
scalapack_wrapper.f90 \
monolis_wrapper_scalapack_c.c

SRC_SOLV_C = \
monolis_solver_c.c \
solver_wrapper.f90

SRC_EIGEN_SOLV_C = \
monolis_eigen_solver_c.c \
eigen_solver_wrapper.f90

SRC_ALL_C = \
$(addprefix define/, $(SRC_DEFINE_C)) \
$(addprefix linalg/, $(SRC_LINALG_C)) \
$(addprefix matrix/, $(SRC_MAT_C)) \
$(addprefix wrapper/, $(SRC_WRAP_C)) \
$(addprefix solver/, $(SRC_SOLV_C)) \
$(addprefix eigen/, $(SRC_EIGEN_SOLV_C))

##> all targes
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

##> lib objs
LIB_SOURCES = \
$(addprefix $(SRC_DIR)/,  $(SRC_ALL)) \
$(addprefix $(WRAP_DIR)/, $(SRC_ALL_C)) \
./src/monolis_solver.f90 \
./src/monolis.f90
LIB_OBJSt   = $(subst $(SRC_DIR), $(OBJ_DIR), $(LIB_SOURCES:.f90=.o))
LIB_OBJS    = $(subst $(WRAP_DIR), $(OBJ_DIR), $(LIB_OBJSt:.c=.o))

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

SRC_SOLVER_C_TEST = \
monolis_solver_c_test.c

SRC_EIGEN_SOLVER_C_TEST = \
monolis_eigen_solver_c_test.c

SRC_ALL_C_TEST = \
$(addprefix define/, $(SRC_DEFINE_C_TEST)) \
$(addprefix linalg/, $(SRC_LINALG_C_TEST)) \
$(addprefix matrix/, $(SRC_MAT_C_TEST)) \
$(addprefix wrapper/, $(SRC_WRAP_C_TEST)) \
$(addprefix solver/, $(SRC_SOLVER_C_TEST)) \
$(addprefix eigen/, $(SRC_EIGEN_SOLVER_C_TEST))

TST_SRC_C_ALL = $(SRC_ALL_C_TEST) monolis_c_test.c
TST_C_SOURCES = $(addprefix $(TST_WRAP_DIR)/, $(TST_SRC_C_ALL))
TST_C_OBJS    = $(subst $(TST_WRAP_DIR), $(OBJ_DIR), $(TST_C_SOURCES:.c=.o))

##> target
all: \
	cp_header \
	cp_header_lib \
	$(LIB_TARGET) \
	$(LIBALL_TARGET) \
	$(TEST_TARGET) \
	$(TEST_C_TARGET)

lib: \
	cp_header \
	cp_header_lib \
	$(LIB_TARGET) \
	$(LIBALL_TARGET)

$(LIB_TARGET): $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS) $(ARC_LIB)

$(TEST_TARGET): $(TST_OBJS)
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) -o $@ $(TST_OBJS) $(USE_LIB)

$(TEST_C_TARGET): $(TST_C_OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(TST_C_OBJS) $(USE_LIB)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(TST_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(WRAP_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(WRAP_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(OBJ_DIR)/%.o: $(TST_WRAP_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

cp_header:
	$(CP) ./wrapper/linalg/monolis_matvec_c.h ./include/
	$(CP) ./wrapper/linalg/monolis_inner_product_c.h ./include/
	$(CP) ./wrapper/linalg/monolis_vec_util_c.h ./include/
	$(CP) ./wrapper/define/monolis_def_struc_c.h ./include/
	$(CP) ./wrapper/define/monolis_def_mat_c.h ./include/
	$(CP) ./wrapper/define/monolis_def_solver_prm_util_c.h ./include/
	$(CP) ./wrapper/define/monolis_def_solver_prm_c.h ./include/
	$(CP) ./wrapper/matrix/monolis_spmat_nzpattern_c.h ./include/
	$(CP) ./wrapper/matrix/monolis_spmat_nzpattern_util_c.h ./include/
	$(CP) ./wrapper/matrix/monolis_spmat_handler_c.h ./include/
	$(CP) ./wrapper/matrix/monolis_spmat_handler_util_c.h ./include/
	$(CP) ./wrapper/matrix/monolis_spmat_copy_c.h ./include/
	$(CP) ./wrapper/wrapper/monolis_wrapper_scalapack_c.h ./include/
	$(CP) ./wrapper/solver/monolis_solver_c.h ./include/
	$(CP) ./wrapper/eigen/monolis_eigen_solver_c.h ./include/
	$(CP) ./wrapper/monolis_solver.h ./include/
	$(CP) ./wrapper/monolis.h ./include/

cp_header_lib:
	$(CP) ./submodule/monolis_utils/bin/* ./bin/
	$(CP) ./submodule/monolis_utils/include/* ./include/
	$(CP) ./submodule/monolis_utils/lib/libmonolis_utils.a ./lib/
	$(CP) ./submodule/gedatsu/bin/* ./bin/
	$(CP) ./submodule/gedatsu/include/* ./include/
	$(CP) ./submodule/gedatsu/lib/libgedatsu.a ./lib/

$(LIBALL_TARGET):
	ar -rc $(LIB_DIR)/libmonolis.a $(LIB_DIR)/libmonolis_solver.a $(LIB_DIR)/libgedatsu.a $(LIB_DIR)/libmonolis_utils.a

clean:
	$(RM) \
	$(LIB_OBJS) \
	$(TST_OBJS) \
	$(TST_C_OBJS) \
	$(LIB_TARGET) \
	$(LIBALL_TARGET) \
	$(TEST_TARGET) \
	$(TEST_C_TARGET) \
	./include/*.mod \
	./include/*.h \
	./bin/*

.PHONY: clean
