#> monolis Makefile

# DIRECTORY setting
INCLUDE  = -I /usr/include -I ./include $(METIS_INC) $(MUMPS_INC)
LIBRARY  = $(METIS_LIB) $(MUMPS_LIB)
BIN_DIR  = ./bin
SRC_DIR  = ./src
OBJ_DIR  = ./obj
LIB_DIR  = ./lib
LIB_LIST = libmonolis.a
MONOLIS_LIB = -L$(LIB_DIR) -lmonolis
CPP      = -cpp $(FLAG_MPI) $(FLAG_METIS) $(FLAG_MUMPS) $(FLAG_TEST) $(FLAG_DEBUG)

# SEC setting
BIN_MS  = gedatsu_simple_mesh_partitioner
#BIN_BC  = gedatsu_bc_partitioner
#BIN_GP  = gedatsu_graph_partitioner
#BIN_CN  = gedatsu_connectivity_partitioner
#BIN_GPR = gedatsu_nodal_value_r_partitioner
#BIN_GPI = gedatsu_nodal_value_i_partitioner
#BIN_CNR = gedatsu_connectivity_value_r_partitioner
#BIN_CNI = gedatsu_connectivity_value_i_partitioner

TARGET_MS = $(addprefix $(BIN_DIR)/, $(BIN_MS))

SRC_MS        = partitioner_simple_mesh.f90

SOURCES_MS = $(addprefix $(SRC_DIR)/, $(SRC_MS))
OBJS_MS = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES_MS:.f90=.o))

# TARGET setting
all: $(TARGET_MS)

$(TARGET_MS): $(OBJS_PART)
	$(FC) $(FFLAGS) -o $@ $(OBJS_PART) $(MONOLIS_LIB) $(LIBRARY)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAG_METIS) -o $@ -c $<

clean:
#	$(RM)

distclean:
#	$(RM)

sampleclean:

.PHONY: clean
