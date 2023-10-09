
CC = mpicc
FC = mpif90
FFLAGS = -O2

TARGET = mesher.out
TARGET2 = a.out

OBJS = mesher.o
OBJS2 = main.o

all: $(TARGET) $(TARGET2)

$(TARGET): $(OBJS)
	$(CC) -o $@ $(OBJS) -L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -lscalapack -lopenblas -L/usr/local/Cellar/gcc/13.2.0/lib/gcc/current/ -lgfortran -L/usr/local/lib -lmpi_mpifh

$(TARGET2): $(OBJS2)
	$(CC) -o $@ $(OBJS2) -L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -lscalapack -lopenblas -L/usr/local/Cellar/gcc/13.2.0/lib/gcc/current/ -lgfortran -L/usr/local/lib -lmpi_mpifh

%.o: %.f90
	$(FC) $(FFLAGS) -I../../include -o $@ -c $<

clean:
	rm mesher.out
	rm a.out
	rm *.o

.PHONY: clean
