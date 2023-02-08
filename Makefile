#> monolis Makefile

FC     = mpif90
FFLAGS = -O2 -mtune=native -march=native -std=legacy -Wno-missing-include-dirs

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
		FFLAGS  = -O2 -std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow -Wno-missing-include-dirs
	endif

	ifeq ($(findstring INTEL, $(DFLAGS)), INTEL)
		FC      = mpiifort
		FFLAGS  = -O2 -align array64byte
		CC      = mpiicc 
		CFLAGS  = -O2 -no-multibyte-chars
		MOD_DIR = -module ./include
	endif

	ifeq ($(findstring TEST, $(DFLAGS)), TEST)
		FLAG_TEST = -DTEST
	endif

	ifeq ($(findstring METIS64, $(DFLAGS)), METIS64)
		FLAG_METIS = -DWITH_METIS64
	endif

	ifeq ($(findstring MUMPS, $(DFLAGS)), MUMPS)
		FLAG_MUMPS = -DWITH_MUMPS
		MUMPS_DIR = .
		MUMPS_INC = -I $(MUMPS_DIR)/include
		MUMPS_LIB = -L$(MUMPS_DIR)/lib -lpord -lmumps_common -ldmumps
	endif

	ifeq ($(findstring OPENMP, $(DFLAGS)), OPENMP)
		FC += -fopenmp
	endif
endif

INCLUDE  = -I /usr/include -I ./include $(METIS_INC) $(MUMPS_INC) -I /Users/morita/git/monolis_utils/include
LIBRARY  = -L/Users/morita/git/monolis_utils/lib -lmonolis_utils $(METIS_LIB) $(MUMPS_LIB)
BIN_DIR  = ./bin
SRC_DIR  = ./src
OBJ_DIR  = ./obj
LIB_DIR  = ./lib
LIB_LIST = libmonolis.a
MONOLIS_LIB = -L$(LIB_DIR) -lmonolis

BIN_PART   = gedatsu_simple_mesh_partitioner
BIN_PARTBC = gedatsu_bc_partitioner
BIN_PARTNG = gedatsu_graph_partitioner
BIN_PARTCN = gedatsu_connectivity_partitioner
BIN_PARTDR = gedatsu_nodal_value_r_partitioner
BIN_PARTDI = gedatsu_nodal_value_i_partitioner
BIN_PARTCR = gedatsu_connectivity_value_r_partitioner
BIN_PARTCI = gedatsu_connectivity_value_i_partitioner
BIN_CONVMG = gedatsu_converter_simple_mesh2graph

CPP      = -cpp $(FLAG_MPI) $(FLAG_METIS) $(FLAG_MUMPS) $(FLAG_TEST) $(FLAG_DEBUG)

MAKE     = make
CD       = cd
RM       = rm -r
AR       = - ar ruv

LIBTARGET  = $(addprefix $(LIB_DIR)/, $(LIB_LIST))
PARTTARGET = $(addprefix $(BIN_DIR)/, $(BIN_PART))
PARTBCTARGET = $(addprefix $(BIN_DIR)/, $(BIN_PARTBC))
PARTNGTARGET = $(addprefix $(BIN_DIR)/, $(BIN_PARTNG))
PARTCNTARGET = $(addprefix $(BIN_DIR)/, $(BIN_PARTCN))
PARTDRTARGET = $(addprefix $(BIN_DIR)/, $(BIN_PARTDR))
PARTDITARGET = $(addprefix $(BIN_DIR)/, $(BIN_PARTDI))
PARTCRTARGET = $(addprefix $(BIN_DIR)/, $(BIN_PARTCR))
PARTCITARGET = $(addprefix $(BIN_DIR)/, $(BIN_PARTCI))
CONVMGTARGET = $(addprefix $(BIN_DIR)/, $(BIN_CONVMG))

SRC_LIST_UTIL   = def_prm.f90 def_mat.f90 def_mesh.f90 util.f90 util_prm.f90
SRC_LIST_MATRIX = fillin.f90 scaling.f90 restruct.f90 matrix_copy.f90 reorder.f90 sparse_util.f90
SRC_LIST_IO     = io_arg.f90 io.f90
SRC_LIST_GRAPH  = graph.f90 graph_comm.f90
SRC_LIST_ALGO   = linalg.f90 matvec.f90 matmat.f90 converge.f90
SRC_LIST_FACT   = 11/fact_LU_11.f90 11/fact_MF_11.f90 33/fact_LU_33.f90 nn/fact_LU_nn.f90 fact_LU.f90 fact_MF.f90
SRC_LIST_PREC   = 33/diag_33.f90 33/sor_33.f90 nn/diag_nn.f90 nn/sor_nn.f90 \
diag.f90 ilu.f90 sor.f90 Jacobi.f90 MUMPS.f90 MF.f90 precond.f90
SRC_LIST_DIRECT = LU.f90
SRC_LIST_ITER   = IR.f90 SOR.f90 CG.f90 GropCG.f90 PipeCR.f90 PipeCG.f90 BiCGSTAB.f90 \
BiCGSTAB_noprec.f90 CABiCGSTAB_noprec.f90 PipeBiCGSTAB.f90 PipeBiCGSTAB_noprec.f90 GMRES.f90
SRC_LIST_SOLVE  = monolis_solve.f90
SRC_LIST_MAIN   = monolis.f90
SRC_LIST_EIGEN  = Lanczos_util.f90 Lanczos.f90 LOBPCG.f90
SRC_LIST_PART   = comm_util.f90 comm_overlap.f90
SRC_LIST_WRAP   = monolis_wrapper_c.c monolis_wrapper.f90

SRC_SOLVER_LIST = \
$(addprefix factorize/, $(SRC_LIST_FACT)) \
$(addprefix precond/, $(SRC_LIST_PREC)) \
$(addprefix direct/, $(SRC_LIST_DIRC)) \
$(addprefix iterative/, $(SRC_LIST_ITER))

SRC_ALL_LIST    = \
$(addprefix util/, $(SRC_LIST_UTIL)) \
$(addprefix io/, $(SRC_LIST_IO)) \
$(addprefix linalg/, $(SRC_LIST_ALGO)) \
$(addprefix graph/, $(SRC_LIST_GRAPH)) \
$(addprefix matrix/, $(SRC_LIST_MATRIX)) \
$(addprefix solver/, $(SRC_SOLVER_LIST)) \
$(addprefix partitioner/, $(SRC_LIST_PART)) \
$(addprefix main/, $(SRC_LIST_SOLVE)) \
$(addprefix eigen/, $(SRC_LIST_EIGEN)) \
$(addprefix wrapper/, $(SRC_LIST_WRAP)) \
$(addprefix main/, $(SRC_LIST_MAIN))

SRC_PART        = partitioner/partitioner_simple_mesh.f90
SRC_PARTBC      = partitioner/partitioner_bc.f90
SRC_PARTNG      = partitioner/partitioner_nodal_graph.f90
SRC_PARTCN      = partitioner/partitioner_connectivity_graph.f90
SRC_PARTDR      = partitioner/partitioner_nodal_val_r.f90
SRC_PARTDI      = partitioner/partitioner_nodal_val_i.f90
SRC_PARTCR      = partitioner/partitioner_connectivity_val_r.f90
SRC_PARTCI      = partitioner/partitioner_connectivity_val_i.f90
SRC_CONVMG      = converter/converter_simple_mesh2graph.f90

SOURCES = $(addprefix $(SRC_DIR)/, $(SRC_ALL_LIST))
OBJSt = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))
OBJS = $(OBJSt:.c=.o)

SOURCES_PART = $(addprefix $(SRC_DIR)/, $(SRC_PART))
OBJS_PART = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PART:.f90=.o))

SOURCES_PARTBC = $(addprefix $(SRC_DIR)/, $(SRC_PARTBC))
OBJS_PARTBC = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PARTBC:.f90=.o))

SOURCES_PARTNG = $(addprefix $(SRC_DIR)/, $(SRC_PARTNG))
OBJS_PARTNG = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PARTNG:.f90=.o))

SOURCES_PARTCN = $(addprefix $(SRC_DIR)/, $(SRC_PARTCN))
OBJS_PARTCN = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PARTCN:.f90=.o))

SOURCES_PARTDR = $(addprefix $(SRC_DIR)/, $(SRC_PARTDR))
OBJS_PARTDR = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PARTDR:.f90=.o))

SOURCES_PARTDI = $(addprefix $(SRC_DIR)/, $(SRC_PARTDI))
OBJS_PARTDI = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PARTDI:.f90=.o))

SOURCES_PARTCR = $(addprefix $(SRC_DIR)/, $(SRC_PARTCR))
OBJS_PARTCR = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PARTCR:.f90=.o))

SOURCES_PARTCI = $(addprefix $(SRC_DIR)/, $(SRC_PARTCI))
OBJS_PARTCI = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_PARTCI:.f90=.o))

SOURCES_CONVMG = $(addprefix $(SRC_DIR)/, $(SRC_CONVMG))
OBJS_CONVMG = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_CONVMG:.f90=.o))

all: $(LIBTARGET) \
	$(PARTTARGET) $(PARTBCTARGET) $(PARTNGTARGET) $(PARTCNTARGET) \
	$(PARTDRTARGET) $(PARTDITARGET) $(PARTCRTARGET) $(PARTCITARGET) \
	$(CONVMGTARGET)

lib: $(LIBTARGET)

$(LIBTARGET): $(OBJS)
	$(AR) $@ $(OBJS)

$(PARTTARGET): $(OBJS_PART)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PART) $(MONOLIS_LIB) $(LIBRARY)

$(PARTBCTARGET): $(OBJS_PARTBC)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PARTBC) $(MONOLIS_LIB) $(LIBRARY)

$(PARTNGTARGET): $(OBJS_PARTNG)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PARTNG) $(MONOLIS_LIB) $(LIBRARY)

$(PARTCNTARGET): $(OBJS_PARTCN)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PARTCN) $(MONOLIS_LIB) $(LIBRARY)

$(PARTDRTARGET): $(OBJS_PARTDR)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PARTDR) $(MONOLIS_LIB) $(LIBRARY)

$(PARTDITARGET): $(OBJS_PARTDI)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PARTDI) $(MONOLIS_LIB) $(LIBRARY)

$(PARTCRTARGET): $(OBJS_PARTCR)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PARTCR) $(MONOLIS_LIB) $(LIBRARY)

$(PARTCITARGET): $(OBJS_PARTCI)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PARTCI) $(MONOLIS_LIB) $(LIBRARY)

$(CONVMGTARGET): $(OBJS_CONVMG)
	$(FC) $(FFLAGS) -o $@ $(OBJS_CONVMG) $(MONOLIS_LIB) $(LIBRARY)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAG_METIS) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(LIBTARGET) \
	$(OBJS_PART) $(OBJS_PARTBC) $(OBJS_PARTNG) $(OBJS_PARTCN) $(OBJS_PARTDR) $(OBJS_PARTDI) $(OBJS_PARTCR) $(OBJS_PARTCI) \
	$(OBJS_CONVMG) \
	./include/*.mod ./bin/*

distclean:
	$(RM) $(OBJS) $(LIBTARGET) \
	/include/*.mod ./bin/*

sampleclean:

.PHONY: clean
