
#!/bin/bash

make clean

make

./mesher.out -i mtx.dat

../../bin/gedatsu_simple_mesh_partitioner -n $1

mpirun -np $1 ./a.out
