
FC = mpif90
FFLAGS = -O2

TARGET = lanczos_driver
OBJS = main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) -L../../lib -lmonolis -lmetis -ldmumps -lmumps_common -lpord -lscalapack -llapack -lblas

%.o: %.f90
	$(FC) $(FFLAGS) -I../../include -o $@ -c $<

clean:
	rm lanczos_driver
	rm *.o

.PHONY: clean
