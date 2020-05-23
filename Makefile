#> monolis Makefile

FC     = mpif90
FFLAGS = -O3 -mtune=native -march=native -mfpmath=both

#CC     = mpicc
#CFLAGS =

ifdef FLAGS
	comma:= ,
	empty:=
	space:= $(empty) $(empty)
	DFLAGS = $(subst $(comma), $(space), $(FLAGS))

	ifeq ($(findstring DEBUG, $(DFLAGS)), DEBUG)
		#FLAG_DEBUG = -DDEBUG
		FFLAGS = -O2 -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow
	endif

	ifeq ($(findstring TEST, $(DFLAGS)), TEST)
		FLAG_TEST = -DTEST
	endif

	ifeq ($(findstring MPI, $(DFLAGS)), MPI)
		FLAG_MPI = -DWITH_MPI
	endif

	ifeq ($(findstring METIS, $(DFLAGS)), METIS)
		FLAG_METIS = -DWITH_METIS
		METIS_DIR = .
		METIS_INC = -I $(METIS_DIR)/include
		METIS_LIB = -L$(METIS_DIR)/lib -lmetis
	endif

	ifeq ($(findstring MUMPS, $(DFLAGS)), MUMPS)
		FLAG_MUMPS = -DWITH_MUMPS
		MUMPS_DIR = .
		MUMPS_INC = -I $(MUMPS_DIR)/include
		MUMPS_LIB = -L$(MUMPS_DIR)/lib -lpord -lmumps_common -ldmumps -lscalapack -lopenblas
	endif
endif

INCLUDE  = -I ./include $(MUMPS_INC)
MOD_DIR  = -J ./include
LIBRARY  = $(METIS_LIB) $(MUMPS_LIB)
BIN_DIR  = ./bin
SRC_DIR  = ./src
OBJ_DIR  = ./obj
LIB_DIR  = ./lib
LIB_LIST = libmonolis.a
MONOLIS_LIB = -L$(LIB_DIR) -lmonolis
BIN_PART = monolis_partitioner
CPP      = -cpp $(FLAG_MPI) $(FLAG_METIS) $(FLAG_MUMPS) $(FLAG_TEST) $(FLAG_DEBUG)

MAKE     = make
CD       = cd
RM       = rm -r
AR       = - ar ruv

LIBTARGET  = $(addprefix $(LIB_DIR)/, $(LIB_LIST))
PARTTARGET = $(addprefix $(BIN_DIR)/, $(BIN_PART))

SRC_LIST_UTIL   = def_prm.f90 def_mat.f90 def_com.f90 def_mesh.f90 util.f90 stdlib.f90 hash.f90
SRC_LIST_MATRIX = fillin.f90 scaling.f90 restruct.f90 reorder.f90 sparse_util.f90
#SRC_LIST_CONV   = convert_full.f90 convert_coo.f90 convert_csr.f90 alloc_matrix.f90
SRC_LIST_IO     = io.f90
SRC_LIST_GRAPH  = graph.f90
SRC_LIST_ALGO   = linalg_com.f90 linalg_util.f90 linalg.f90 matvec.f90 matmat.f90 converge.f90
SRC_LIST_FACT   = 11/fact_LU_11.f90 33/fact_LU_33.f90 nn/fact_LU_nn.f90 fact_LU.f90
SRC_LIST_PREC   = 33/diag_33.f90 33/sor_33.f90 nn/diag_nn.f90 nn/sor_nn.f90 diag.f90 ilu.f90 sor.f90 Jacobi.f90 MUMPS.f90 precond.f90
SRC_LIST_DIRECT = LU.f90
SRC_LIST_ITER   = IR.f90 SOR.f90 CG.f90 GropCG.f90 PipeCR.f90 PipeCG.f90 BiCGSTAB.f90 BiCGSTAB_noprec.f90 CABiCGSTAB_noprec.f90 PipeBiCGSTAB.f90 PipeBiCGSTAB_noprec.f90
SRC_LIST_MAIN   = monolis_solve.f90 monolis.f90
SRC_LIST_PART   = comm_util.f90 comm_overlap.f90
SRC_LIST_PART_BIN = partitioner.f90
#SRC_LIST_WRAP  =

SRC_SOLVER_LIST = $(addprefix factorize/, $(SRC_LIST_FACT)) $(addprefix precond/, $(SRC_LIST_PREC)) $(addprefix direct/, $(SRC_LIST_DIRC)) $(addprefix iterative/, $(SRC_LIST_ITER))
SRC_ALL_LIST    = $(addprefix util/, $(SRC_LIST_UTIL)) $(addprefix io/, $(SRC_LIST_IO)) $(addprefix graph/, $(SRC_LIST_GRAPH))  $(addprefix linalg/, $(SRC_LIST_ALGO)) $(addprefix matrix/, $(SRC_LIST_MATRIX)) $(addprefix solver/, $(SRC_SOLVER_LIST)) $(addprefix partitioner/, $(SRC_LIST_PART)) $(addprefix main/, $(SRC_LIST_MAIN))
SRC_PART_LIST   = $(addprefix partitioner/, $(SRC_LIST_PART_BIN)) $(SRC_ALL_LIST)

SOURCES = $(addprefix $(SRC_DIR)/, $(SRC_ALL_LIST))
OBJS = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))

SOURCES_PART = $(addprefix $(SRC_DIR)/, $(SRC_PART_LIST))
OBJS_PART = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PART:.f90=.o))

all: $(LIBTARGET) $(PARTTARGET)
#	$(CD) sample && $(MAKE)

$(LIBTARGET): $(OBJS)
	$(AR) $@ $(OBJS)

$(PARTTARGET): $(OBJS_PART)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PART) $(LIBRARY)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(LIBTARGET) $(PARTTARGET) ./include/*.mod

distclean:
	$(RM) $(OBJS) $(LIBTARGET) $(PARTTARGET) ./include/*.mod

sampleclean:

.PHONY: clean
