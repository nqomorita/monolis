
FC       = gfortran
FFLAGS   = -O2 -fbounds-check -fbacktrace -ffpe-trap=invalid
LDFLAGS  =
LIBS     =
INCLUDE  = -I ./include
MOD_DIR  = -J ./include
BIN_DIR  = ./bin
SRC_DIR  = ./src
OBJ_DIR  = ./obj
LIB_DIR  = ./lib
BIN_LIST = monolis
LIB_LIST = libmonolis.a
TARGET   = $(addprefix $(BIN_DIR)/, $(BIN_LIST))
LIBTARGET= $(addprefix $(LIB_DIR)/, $(LIB_LIST))
SRC_LIST = def_prm.f90 def_mat.f90 def_com.f90 util.f90 hecmw_monolis.f90 monolis.f90 main.f90
LIB_SRC_LIST = $(subst main.f90, , $(SRC_LIST))
SOURCES  = $(addprefix $(SRC_DIR)/, $(SRC_LIST))
LIB_SOURCES  = $(addprefix $(SRC_DIR)/, $(LIB_SRC_LIST))
OBJS     = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))
LIB_OBJS = $(subst $(SRC_DIR), $(OBJ_DIR), $(LIB_SOURCES:.f90=.o))
RM       = rm
AR       = - ar ruv

all: $(TARGET) $(LIBTARGET)

$(TARGET): $(OBJS) $(LIBS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

$(LIBTARGET): $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET) ./include/*.mod

.PHONY: clean
