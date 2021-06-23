#> monolis Makefile

FC     = mpif90
FFLAGS = -O2 -mtune=native -march=native -std=legacy

CC     = mpic++
CFLAGS = -O2

MOD_DIR  = -J ./include

# DEFAULT FLAGS
FLAG_MPI   = -DWITH_MPI
FLAG_METIS = -DWITH_METIS
METIS_DIR = .
METIS_INC = -I $(METIS_DIR)/include
METIS_LIB = -L$(METIS_DIR)/lib -lmetis

ifdef FLAGS
	comma:= ,
	empty:=
	space:= $(empty) $(empty)
	DFLAGS = $(subst $(comma), $(space), $(FLAGS))

	ifeq ($(findstring DEBUG, $(DFLAGS)), DEBUG)
		FFLAGS  = -O2 -std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow
	endif

	ifeq ($(findstring INTEL, $(DFLAGS)), INTEL)
		FC      = mpiifort
		FFLAGS  = -O2 -align array64byte
		CC      = mpiicpc
		CFLAGS  = -O2
		MOD_DIR = -module ./include
	endif

	ifeq ($(findstring TEST, $(DFLAGS)), TEST)
		FLAG_TEST = -DTEST
	endif

#	ifeq ($(findstring MPI, $(DFLAGS)), MPI)
#		FLAG_MPI = -DWITH_MPI
#	endif

#	ifeq ($(findstring METIS, $(DFLAGS)), METIS)
#		FLAG_METIS = -DWITH_METIS
#		METIS_DIR = .
#		METIS_INC = -I $(METIS_DIR)/include
#		METIS_LIB = -L$(METIS_DIR)/lib -lmetis
#	endif

	ifeq ($(findstring MUMPS, $(DFLAGS)), MUMPS)
		FLAG_MUMPS = -DWITH_MUMPS
		MUMPS_DIR = .
		MUMPS_INC = -I $(MUMPS_DIR)/include
		MUMPS_LIB = -L$(MUMPS_DIR)/lib -lpord -lmumps_common -ldmumps
	endif
endif

INCLUDE  = -I ./include $(METIS_INC) -I ./include $(MUMPS_INC)
LIBRARY  = $(METIS_LIB) $(MUMPS_LIB)
BIN_DIR  = ./bin
SRC_DIR  = ./src
OBJ_DIR  = ./obj
LIB_DIR  = ./lib
LIB_LIST = libmonolis.a
MONOLIS_LIB = -L$(LIB_DIR) -lmonolis
BIN_PART = monolis_partitioner
BIN_PARTBC = monolis_bc_partitioner
BIN_REF1 = monolis_h_refiner_hex
BIN_REF2 = monolis_h_refiner_tet
BIN_REF3 = monolis_p_refiner_hex
BIN_REF4 = monolis_p_refiner_tet
BIN_DBC1 = monolis_dbc_all_surf_tet
BIN_DBC2 = monolis_dbc_all_surf_hex
BIN_TEST = monolis_test
CPP      = -cpp $(FLAG_MPI) $(FLAG_METIS) $(FLAG_MUMPS) $(FLAG_TEST) $(FLAG_DEBUG)

MAKE     = make
CD       = cd
RM       = rm -r
AR       = - ar ruv

LIBTARGET  = $(addprefix $(LIB_DIR)/, $(LIB_LIST))
PARTTARGET = $(addprefix $(BIN_DIR)/, $(BIN_PART))
PARTBCTARGET = $(addprefix $(BIN_DIR)/, $(BIN_PARTBC))
REF1TARGET = $(addprefix $(BIN_DIR)/, $(BIN_REF1))
REF2TARGET = $(addprefix $(BIN_DIR)/, $(BIN_REF2))
REF3TARGET = $(addprefix $(BIN_DIR)/, $(BIN_REF3))
REF4TARGET = $(addprefix $(BIN_DIR)/, $(BIN_REF4))
DBC1TARGET = $(addprefix $(BIN_DIR)/, $(BIN_DBC1))
DBC2TARGET = $(addprefix $(BIN_DIR)/, $(BIN_DBC2))
TESTTARGET = $(addprefix $(BIN_DIR)/, $(BIN_TEST))

SRC_LIST_UTIL   = def_prm.f90 def_mat.f90 def_com.f90 def_mesh.f90 util.f90 stdlib.f90 hash.f90
SRC_LIST_MATRIX = fillin.f90 scaling.f90 restruct.f90 matrix_copy.f90 reorder.f90 sparse_util.f90
#SRC_LIST_CONV   = convert_full.f90 convert_coo.f90 convert_csr.f90 alloc_matrix.f90
SRC_LIST_IO     = io.f90 io_mtx.f90
SRC_LIST_GRAPH  = graph.f90
SRC_LIST_SHAPE  = shape_util.f90 shape_C2D3.f90 shape_C2D4.f90 shape_C2D6.f90 shape_C3D4.f90 shape_C3D8.f90 shape_C3D10.f90
SRC_LIST_GEOM   = geom.f90 neighbor_search.f90
SRC_LIST_ALGO   = linalg_com.f90 linalg_util.f90 linalg.f90 matvec.f90 matmat.f90 converge.f90
SRC_LIST_FACT   = 11/fact_LU_11.f90 11/fact_MF_11.f90 33/fact_LU_33.f90 nn/fact_LU_nn.f90 fact_LU.f90 fact_MF.f90
SRC_LIST_PREC   = 33/diag_33.f90 33/sor_33.f90 33/ROM_GS_33.f90 nn/diag_nn.f90 nn/sor_nn.f90 diag.f90 ilu.f90 sor.f90 Jacobi.f90 MUMPS.f90 ROM_GS.f90 ROM.f90 MF.f90 precond.f90
SRC_LIST_DIRECT = LU.f90
SRC_LIST_ITER   = IR.f90 SOR.f90 CG.f90 GropCG.f90 PipeCR.f90 PipeCG.f90 BiCGSTAB.f90 BiCGSTAB_noprec.f90 CABiCGSTAB_noprec.f90 PipeBiCGSTAB.f90 PipeBiCGSTAB_noprec.f90 GMRES.f90
SRC_LIST_SOLVE  = monolis_solve.f90
SRC_LIST_MAIN   = monolis.f90
SRC_LIST_EIGEN  = Lanczos_util.f90 Lanczos.f90 LOBPCG.f90
SRC_LIST_PART   = comm_util.f90 comm_overlap.f90
SRC_LIST_REFN   = p_refiner.f90 dbc_all_util.f90
SRC_LIST_WRAP   = monolis_wrapper_c.c monolis_wrapper.f90

SRC_SOLVER_LIST = $(addprefix factorize/, $(SRC_LIST_FACT)) $(addprefix precond/, $(SRC_LIST_PREC)) $(addprefix direct/, $(SRC_LIST_DIRC)) $(addprefix iterative/, $(SRC_LIST_ITER))
SRC_ALL_LIST    = $(addprefix util/, $(SRC_LIST_UTIL)) $(addprefix io/, $(SRC_LIST_IO)) $(addprefix graph/, $(SRC_LIST_GRAPH)) $(addprefix shape/, $(SRC_LIST_SHAPE)) $(addprefix geom/, $(SRC_LIST_GEOM)) $(addprefix linalg/, $(SRC_LIST_ALGO)) $(addprefix matrix/, $(SRC_LIST_MATRIX)) $(addprefix solver/, $(SRC_SOLVER_LIST)) $(addprefix partitioner/, $(SRC_LIST_PART)) $(addprefix refiner/, $(SRC_LIST_REFN)) $(addprefix main/, $(SRC_LIST_SOLVE)) $(addprefix eigen/, $(SRC_LIST_EIGEN)) $(addprefix wrapper/, $(SRC_LIST_WRAP)) $(addprefix main/, $(SRC_LIST_MAIN))
SRC_PART        = partitioner/partitioner.f90
SRC_PARTBC      = partitioner/partitioner_bc.f90
SRC_REF1        = refiner/h_refiner_hex.f90
SRC_REF2        = refiner/h_refiner_tet.f90
SRC_REF3        = refiner/p_refiner_hex.f90
SRC_REF4        = refiner/p_refiner_tet.f90
SRC_DBC1        = refiner/dbc_all_tet.f90
SRC_DBC2        = refiner/dbc_all_hex.f90
SRC_TEST        = util/test.f90

SOURCES = $(addprefix $(SRC_DIR)/, $(SRC_ALL_LIST))
OBJSt = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))
OBJS = $(OBJSt:.c=.o)

SOURCES_PART = $(addprefix $(SRC_DIR)/, $(SRC_PART))
OBJS_PART = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PART:.f90=.o))

SOURCES_PARTBC = $(addprefix $(SRC_DIR)/, $(SRC_PARTBC))
OBJS_PARTBC = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PARTBC:.f90=.o))

SOURCES_REF1 = $(addprefix $(SRC_DIR)/, $(SRC_REF1))
OBJS_REF1 = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_REF1:.f90=.o))

SOURCES_REF2 = $(addprefix $(SRC_DIR)/, $(SRC_REF2))
OBJS_REF2 = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_REF2:.f90=.o))

SOURCES_REF3 = $(addprefix $(SRC_DIR)/, $(SRC_REF3))
OBJS_REF3 = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_REF3:.f90=.o))

SOURCES_REF4 = $(addprefix $(SRC_DIR)/, $(SRC_REF4))
OBJS_REF4 = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_REF4:.f90=.o))

SOURCES_DBC1 = $(addprefix $(SRC_DIR)/, $(SRC_DBC1))
OBJS_DBC1 = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_DBC1:.f90=.o))

SOURCES_DBC2 = $(addprefix $(SRC_DIR)/, $(SRC_DBC2))
OBJS_DBC2 = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_DBC2:.f90=.o))

SOURCES_TEST = $(addprefix $(SRC_DIR)/, $(SRC_TEST))
OBJS_TEST = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_TEST:.f90=.o))

all: $(LIBTARGET) $(PARTTARGET) $(PARTBCTARGET) $(REF1TARGET) $(REF2TARGET) $(REF3TARGET) $(REF4TARGET) \
$(DBC1TARGET) $(DBC2TARGET) $(TESTTARGET)

$(LIBTARGET): $(OBJS)
	$(AR) $@ $(OBJS)

$(PARTTARGET): $(OBJS_PART)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PART) $(MONOLIS_LIB) $(LIBRARY)

$(PARTBCTARGET): $(OBJS_PARTBC)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PARTBC) $(MONOLIS_LIB) $(LIBRARY)

$(REF1TARGET): $(OBJS_REF1)
	$(FC) $(FFLAGS) -o $@ $(OBJS_REF1) $(MONOLIS_LIB) $(LIBRARY)

$(REF2TARGET): $(OBJS_REF2)
	$(FC) $(FFLAGS) -o $@ $(OBJS_REF2) $(MONOLIS_LIB) $(LIBRARY)

$(REF3TARGET): $(OBJS_REF3)
	$(FC) $(FFLAGS) -o $@ $(OBJS_REF3) $(MONOLIS_LIB) $(LIBRARY)

$(REF4TARGET): $(OBJS_REF4)
	$(FC) $(FFLAGS) -o $@ $(OBJS_REF4) $(MONOLIS_LIB) $(LIBRARY)

$(DBC1TARGET): $(OBJS_DBC1)
	$(FC) $(FFLAGS) -o $@ $(OBJS_DBC1) $(MONOLIS_LIB) $(LIBRARY)

$(DBC2TARGET): $(OBJS_DBC2)
	$(FC) $(FFLAGS) -o $@ $(OBJS_DBC2) $(MONOLIS_LIB) $(LIBRARY)

$(TESTTARGET): $(OBJS_TEST)
	$(FC) $(FFLAGS) -o $@ $(OBJS_TEST) $(MONOLIS_LIB) $(LIBRARY)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAG_METIS) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(LIBTARGET) $(PARTTARGET) $(REF1TARGET) $(REF2TARGET) $(REF3TARGET) $(REF4TARGET) \
	$(DBC1TARGET) $(DBC2TARGET) $(TESTTARGET) ./include/*.mod

distclean:
	$(RM) $(OBJS) $(LIBTARGET) $(PARTTARGET) $(REF1TARGET) $(REF2TARGET) $(REF3TARGET) $(REF4TARGET) \
	$(DBC1TARGET) $(DBC2TARGET) $(TESTTARGET) /include/*.mod

sampleclean:

.PHONY: clean
