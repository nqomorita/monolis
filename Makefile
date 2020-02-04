
FLAG_MPI   = -DWITH_MPI
FLAG_METIS = -DWITH_METIS
FLAG_MUMPS = -DWITH_MUMPS
FLAG_DDM   = -DOVER_DDM
#FLAG_DEBUG = -DDEBUG
#FLAG_TEST  = -DTEST_ALL
export CPP        = -cpp $(FLAG_MPI) $(FLAG_METIS) $(FLAG_MUMPS) $(FLAG_TEST) $(FLAG_DEBUG)

export FC         = mpif90
export FFLAGS     = -O2 -fbounds-check -fbacktrace -ffpe-trap=invalid
export CC         = mpicc
export CFLAGS     =

export METIS_DIR  = $(HOME)/FrontISTR_tools
export METIS_INC  = -I $(METIS_DIR)/include
export METIS_LIB  = -L$(METIS_DIR)/lib -lmetis

export MUMPS_DIR  = $(HOME)/FrontISTR_tools
export MUMPS_INC  = -I $(MUMPS_DIR)/include
export MUMPS_LIB  = -L$(MUMPS_DIR)/lib -lpord -lmumps_common -ldmumps -lscalapack -lopenblas

export INCLUDE    = -I ./include $(MUMPS_INC)
export MOD_DIR    = -J ./include
export LIBRARY    = $(METIS_LIB) $(MUMPS_LIB)
export SRC_DIR    = ./src
export SMP_DIR    = ./sample
export OBJ_DIR    = ./obj
export LIB_DIR    = ./lib
export LIB_LIST   = libmonolis.a

export MAKE       = make
export CD         = cd
export RM         = rm -r
export AR         = - ar ruv

LIBTARGET  = $(addprefix $(LIB_DIR)/, $(LIB_LIST))

SRC_LIST_UTIL = def_prm.f90 def_mat.f90 def_com.f90 util.f90 fillin.f90 hash.f90
#SRC_LIST_CONV = convert_full.f90 convert_coo.f90 convert_csr.f90 alloc_matrix.f90
SRC_LIST_ALGO = linalg_com.f90 linalg_util.f90 linalg.f90 matvec.f90 matmat.f90 converge.f90 scaling.f90 restruct.f90 reorder.f90
SRC_LIST_FACT = 11/fact_LU_11.f90 33/fact_LU_33.f90 nn/fact_LU_nn.f90 fact_LU.f90
SRC_LIST_PREC = 33/diag_33.f90 33/sor_33.f90 nn/diag_nn.f90 nn/sor_nn.f90 diag.f90 ilu.f90 sor.f90 Jacobi.f90 MUMPS.f90 precond.f90
SRC_LIST_DIRC = LU.f90
SRC_LIST_ITER = IR.f90 SOR.f90 CG.f90 GropCG.f90 PipeCR.f90 PipeCG.f90 BiCGSTAB.f90 BiCGSTAB_noprec.f90 CABiCGSTAB_noprec.f90 PipeBiCGSTAB.f90 PipeBiCGSTAB_noprec.f90
SRC_LIST_LIB  = monolis_solve.f90 monolis.f90

SRC_ALL_LIST  = $(addprefix util/, $(SRC_LIST_UTIL)) $(addprefix convert/, $(SRC_LIST_CONV)) $(addprefix linalg/, $(SRC_LIST_ALGO)) $(addprefix factorize/, $(SRC_LIST_FACT)) $(addprefix precond/, $(SRC_LIST_PREC)) $(addprefix direct/, $(SRC_LIST_DIRC)) $(addprefix iterative/, $(SRC_LIST_ITER)) $(addprefix main/, $(SRC_LIST_LIB))

SOURCES = $(addprefix $(SRC_DIR)/, $(SRC_ALL_LIST))

OBJS = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))

all: $(LIBTARGET)
	$(CD) sample && $(MAKE)

$(LIBTARGET): $(OBJS)
	$(AR) $@ $(OBJS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(LIBTARGET) ./include/*.mod

distclean:
	$(RM) $(OBJS) $(LIBTARGET) ./include/*.mod

sampleclean:

.PHONY: clean
