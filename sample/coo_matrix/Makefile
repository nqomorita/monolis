
TARGET = monolis_sample
OBJS = main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBRARY) -L../../lib -lmonolis

%.o: %.f90
	$(FC) $(FFLAGS) $(CPP) -I../../include -o $@ -c $<

clean:

distclean:

sampleclean:

.PHONY: clean
