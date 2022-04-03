
TARGET = main

OBJS   = main.o
FC     = mpif90
CC     = mpicc
CPP    = -cpp

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) -L../../lib -lmonolis -lmetis -llapack

%.o: %.f90
	$(FC) $(FFLAGS) $(CPP) -I../../include -o $@ -c $<

%.o: %.c
	$(CC) $(CFLAGS) -I../../include -o $@ -c $<

clean:
	rm main.o main

distclean:

sampleclean:

.PHONY: clean
