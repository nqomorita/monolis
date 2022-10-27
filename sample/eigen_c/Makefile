
FC = mpif90
FFLAGS = -O2

CC = mpic++
CFLAGS = -O2

TARGET = eigen_test
OBJS = main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) -L../../lib -lmonolis -lmetis -llapack

%.o: %.c
	$(CC) $(CFLAGS) -I../../include -o $@ -c $<

clean:
	rm eigen_test
	rm *.o

.PHONY: clean
