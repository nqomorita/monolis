
FC = mpif90
FFLAGS = -O2

CC = mpicc
CFLAGS = -O2

TARGET = solver_test
OBJS = main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) -I../../include -o solver_test $(OBJS) -L../../lib -lmonolis -lmetis -llapack

%.o: %.c
	$(CC) $(CFLAGS) -I../../include -o $@ -c $<

clean:
	rm solver_test; rm *.o

.PHONY: clean
