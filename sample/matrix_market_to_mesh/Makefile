
TARGET1 = matrix_market_to_mesh
TARGET2 = matrix_market_solve

OBJS1  = main1.o
OBJS2  = main2.o
FC     = mpif90
CPP    = -cpp

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJS1)
	$(FC) -o $@ $(OBJS1) -L../../lib -lmonolis -lmetis -llapack

$(TARGET2): $(OBJS2)
	$(FC) -o $@ $(OBJS2) -L../../lib -lmonolis -lmetis -llapack

%.o: %.f90
	$(FC) $(FFLAGS) $(CPP) -I../../include -o $@ -c $<

clean:
	rm main1.o main2.o matrix_market_to_mesh matrix_market_solve

distclean:

sampleclean:

.PHONY: clean
