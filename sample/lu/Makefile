FC = mpif90
FFLAGS = -O2

TARGET = mesher.out
TARGET2 = a.out
OBJS = mesher.o
OBJS2 = main.o

all: $(TARGET) $(TARGET2)

LIBS = -L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -ldmumps -lmumps_common -lpord -lscalapack -llapack -lblas

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS)  $(LIBS)

$(TARGET2): $(OBJS2)
	$(FC) -o $@ $(OBJS2) $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -I../../include -o $@ -c $<

clean:
	rm mesher.out
	rm a.out
	rm *.o

.PHONY: clean

