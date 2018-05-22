
FLAG_MPI  = -DWITHMPI
CPP       = -cpp $(FLAG_MPI)

FC        = mpif90
FFLAGS    = -O2 -fbounds-check -fbacktrace -ffpe-trap=invalid
LDFLAGS   =
LIBS      =

INCLUDE   = -I ./include
MOD_DIR   = -J ./include
BIN_DIR   = ./bin
SRC_DIR   = ./src
OBJ_DIR   = ./obj
LIB_DIR   = ./lib
BIN_LIST  = monolis
LIB_LIST  = libmonolis.a
RM        = rm -r
AR        = - ar ruv

TARGET    = $(addprefix $(BIN_DIR)/, $(BIN_LIST))
LIBTARGET = $(addprefix $(LIB_DIR)/, $(LIB_LIST))

SRC_LIST_UTIL = def_prm.f90 def_mat.f90 def_com.f90 util.f90 hecmw_monolis.f90
SRC_LIST_ALGO = linalg_com.f90 linalg.f90
SRC_LIST_PREC = precond.f90
SRC_LIST_ITER = IR.f90 iterative.f90
SRC_LIST_LIB  = monolis.f90
SRC_LIST_MAIN = main.f90

SRC_LIST    = $(addprefix util/, $(SRC_LIST_UTIL)) $(addprefix linalg/, $(SRC_LIST_ALGO)) $(addprefix precond/, $(SRC_LIST_PREC)) $(addprefix iterative/, $(SRC_LIST_ITER)) $(addprefix main/, $(SRC_LIST_LIB)) $(addprefix main/, $(SRC_LIST_MAIN))
SRC_LIST_AR = $(addprefix util/, $(SRC_LIST_UTIL)) $(addprefix linalg/, $(SRC_LIST_ALGO)) $(addprefix precond/, $(SRC_LIST_PREC)) $(addprefix iterative/, $(SRC_LIST_ITER)) $(addprefix main/, $(SRC_LIST_LIB))

SOURCES    = $(addprefix $(SRC_DIR)/, $(SRC_LIST))
SOURCES_AR = $(addprefix $(SRC_DIR)/, $(SRC_LIST_AR))

OBJS    = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))
OBJS_AR = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_AR:.f90=.o))

all: $(TARGET) $(LIBTARGET)

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

$(LIBTARGET): $(OBJS_AR)
	$(AR) $@ $(OBJS_AR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET) $(LIBTARGET) ./include/*.mod

.PHONY: clean
