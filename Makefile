
FLAG_MPI   = -DWITH_MPI
FLAG_METIS = -DWITH_METIS
FLAG_MUMPS = -DWITH_MUMPS
FLAG_DDM   = -DOVER_DDM
#FLAG_DEBUG = -DDEBUG
#FLAG_TEST  = -DTEST_ALL
CPP        = -cpp $(FLAG_MPI) $(FLAG_METIS) $(FLAG_MUMPS) $(FLAG_TEST) $(FLAG_DEBUG)

FC         = mpif90
FFLAGS     = -O2 -fbounds-check -fbacktrace -ffpe-trap=invalid
CC         = mpicc
CFLAGS     =

METIS_DIR  = $(HOME)/FrontISTR_tools
METIS_INC  = -I $(METIS_DIR)/include
METIS_LIB  = -L$(METIS_DIR)/lib -lmetis

MUMPS_DIR  = $(HOME)/FrontISTR_tools
MUMPS_INC  = -I $(MUMPS_DIR)/include
MUMPS_LIB  = -L$(MUMPS_DIR)/lib -lpord -lmumps_common -ldmumps -lscalapack -lopenblas

INCLUDE    = -I ./include $(MUMPS_INC)
MOD_DIR    = -J ./include
LIBRARY    = $(METIS_LIB) $(MUMPS_LIB)
BIN_DIR    = ./bin
SRC_DIR    = ./src
SMP_DIR    = ./sample
OBJ_DIR    = ./obj
LIB_DIR    = ./lib
BIN_LIST   = monolis
LIB_LIST   = libmonolis.a
SMP1_LIST  = hash_table/monolis_sample
SMP2_LIST  = coo_matrix/monolis_sample
SMP3_LIST  = coo_matrix_c/monolis_sample
SMP4_LIST  = csr_matrix/monolis_sample
SMP5_LIST  = csr_matrix_c/monolis_sample
RM         = rm -r
AR         = - ar ruv

TARGET     = $(addprefix $(BIN_DIR)/, $(BIN_LIST))
LIBTARGET  = $(addprefix $(LIB_DIR)/, $(LIB_LIST))
SMP1TARGET  = $(addprefix $(SMP_DIR)/, $(SMP1_LIST))
SMP2TARGET  = $(addprefix $(SMP_DIR)/, $(SMP2_LIST))
SMP3TARGET  = $(addprefix $(SMP_DIR)/, $(SMP3_LIST))
SMP4TARGET  = $(addprefix $(SMP_DIR)/, $(SMP4_LIST))
SMP5TARGET  = $(addprefix $(SMP_DIR)/, $(SMP5_LIST))

SRC_LIST_UTIL = def_prm.f90 def_mat.f90 def_com.f90 util.f90 fillin.f90 hash.f90
#SRC_LIST_CONV = convert_full.f90 convert_coo.f90 convert_csr.f90 alloc_matrix.f90
SRC_LIST_ALGO = linalg_com.f90 linalg_util.f90 linalg.f90 matvec.f90 matmat.f90 converge.f90 scaling.f90 restruct.f90 reorder.f90
SRC_LIST_FACT = 11/fact_LU_11.f90 33/fact_LU_33.f90 nn/fact_LU_nn.f90 fact_LU.f90
SRC_LIST_PREC = 33/diag_33.f90 33/sor_33.f90 nn/diag_nn.f90 nn/sor_nn.f90 diag.f90 ilu.f90 sor.f90 Jacobi.f90 MUMPS.f90 precond.f90
SRC_LIST_DIRC = LU.f90
SRC_LIST_ITER = IR.f90 SOR.f90 CG.f90 GropCG.f90 PipeCR.f90 PipeCG.f90 BiCGSTAB.f90 BiCGSTAB_noprec.f90 CABiCGSTAB_noprec.f90 PipeBiCGSTAB.f90 PipeBiCGSTAB_noprec.f90
SRC_LIST_LIB  = monolis_solve.f90 monolis.f90
SRC_LIST_MAIN = main.f90
SRC_LIST_SMP1 = hash_table/main.f90
SRC_LIST_SMP2 = coo_matrix/main.f90
SRC_LIST_SMP3 = coo_matrix_c/main.c
SRC_LIST_SMP4 = csr_matrix/main.f90
SRC_LIST_SMP5 = csr_matrix_c/main.c

SRC_ALL_LIST    = $(addprefix util/, $(SRC_LIST_UTIL)) $(addprefix convert/, $(SRC_LIST_CONV)) $(addprefix linalg/, $(SRC_LIST_ALGO)) $(addprefix factorize/, $(SRC_LIST_FACT)) $(addprefix precond/, $(SRC_LIST_PREC)) $(addprefix direct/, $(SRC_LIST_DIRC)) $(addprefix iterative/, $(SRC_LIST_ITER)) $(addprefix main/, $(SRC_LIST_LIB)) $(addprefix main/, $(SRC_LIST_MAIN))
SRC_ALL_LIST_AR = $(addprefix util/, $(SRC_LIST_UTIL)) $(addprefix convert/, $(SRC_LIST_CONV)) $(addprefix linalg/, $(SRC_LIST_ALGO)) $(addprefix factorize/, $(SRC_LIST_FACT)) $(addprefix precond/, $(SRC_LIST_PREC)) $(addprefix direct/, $(SRC_LIST_DIRC)) $(addprefix iterative/, $(SRC_LIST_ITER)) $(addprefix main/, $(SRC_LIST_LIB))

SOURCES    = $(addprefix $(SRC_DIR)/, $(SRC_ALL_LIST))
SOURCES_AR = $(addprefix $(SRC_DIR)/, $(SRC_ALL_LIST_AR))
SAMPLE1    = $(addprefix $(SMP_DIR)/, $(SRC_LIST_SMP1))
SAMPLE2    = $(addprefix $(SMP_DIR)/, $(SRC_LIST_SMP2))
SAMPLE3    = $(addprefix $(SMP_DIR)/, $(SRC_LIST_SMP3))
SAMPLE4    = $(addprefix $(SMP_DIR)/, $(SRC_LIST_SMP4))
SAMPLE5    = $(addprefix $(SMP_DIR)/, $(SRC_LIST_SMP5))

OBJS    = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))
OBJS_AR = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_AR:.f90=.o))
SMP1    = $(SAMPLE1:.f90=.o)
SMP2    = $(SAMPLE2:.f90=.o)
SMP3    = $(SAMPLE3:.c=.o)
SMP4    = $(SAMPLE4:.f90=.o)
SMP5    = $(SAMPLE5:.c=.o)

all: $(TARGET) $(LIBTARGET) $(SMP1TARGET) $(SMP2TARGET) $(SMP3TARGET) $(SMP4TARGET) $(SMP5TARGET)

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBRARY)

$(LIBTARGET): $(OBJS_AR)
	$(AR) $@ $(OBJS_AR)

$(SMP1TARGET): $(SMP1)
	$(FC) -o $@ $(SMP1) $(LIBRARY) -L$(LIB_DIR) -lmonolis

$(SMP2TARGET): $(SMP2)
	$(FC) -o $@ $(SMP2) $(LIBRARY) -L$(LIB_DIR) -lmonolis

$(SMP3TARGET): $(SMP3)
	$(FC) -o $@ $(SMP3) $(LIBRARY) -L$(LIB_DIR) -lmonolis

$(SMP4TARGET): $(SMP4)
	$(FC) -o $@ $(SMP4) $(LIBRARY) -L$(LIB_DIR) -lmonolis

$(SMP5TARGET): $(SMP5)
	$(FC) -o $@ $(SMP5) $(LIBRARY) -L$(LIB_DIR) -lmonolis

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(SMP_DIR)/%.o: $(SMP_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) -o $@ -c $<

$(SMP_DIR)/%.o: $(SMP_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(SMP1) $(SMP2) $(SMP3) $(SMP4) $(SMP5) $(TARGET) $(LIBTARGET) $(SMP1TARGET) $(SMP2TARGET) $(SMP3TARGET) $(SMP4TARGET) $(SMP5TARGET) ./include/*.mod

distclean:
	$(RM) $(OBJS) $(SMP1) $(SMP2) $(SMP3) $(SMP4) $(SMP5) $(TARGET) $(LIBTARGET) $(SMP1TARGET) $(SMP2TARGET) $(SMP3TARGET) $(SMP4TARGET) $(SMP5TARGET) ./include/*.mod

sampleclean:
	$(RM) $(SMP1) $(SMP2) $(SMP3) $(SMP4) $(SMP5) $(SMP1TARGET) $(SMP2TARGET) $(SMP3TARGET) $(SMP4TARGET) $(SMP5TARGET)

.PHONY: clean
