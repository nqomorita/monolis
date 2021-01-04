
TARGET = matrix_market_to_mesh
OBJS   = main.o
FC     = mpif90
CPP    = -cpp

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) -L../../lib -lmonolis -lmetis -llapack

%.o: %.f90
	$(FC) $(FFLAGS) $(CPP) -I../../include -o $@ -c $<

clean:
	rm main.o

distclean:

sampleclean:

.PHONY: clean
