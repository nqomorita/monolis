
TARGET = lobpcg
OBJS = lobpcg.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBRARY) -L../../lib -lmonolis -lmetis -llapack

%.o: %.f90
	$(FC) $(FFLAGS) $(CPP) -I../../include -o $@ -c $<

clean:

distclean:

sampleclean:

.PHONY: clean
