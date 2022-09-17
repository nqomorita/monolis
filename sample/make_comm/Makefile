
TARGET1 = make_comm

OBJS1  = main.o
FC     = mpif90 -std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow
CPP    = -cpp

all: $(TARGET1)

$(TARGET1): $(OBJS1)
	$(FC) -o $@ $(OBJS1) -L../../lib -lmonolis -lmetis

%.o: %.f90
	$(FC) $(FFLAGS) $(CPP) -I../../include -o $@ -c $<

clean:
	rm main.o make_comm

distclean:

sampleclean:

.PHONY: clean
