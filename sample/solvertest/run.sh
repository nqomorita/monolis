
#!/bin/bash

echo "./a.out -i young1c.mtx"

mpif90 -I../../include -std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow mesher.f90 -L../../lib -lmonolis -lmetis

./a.out -i bcsstk14.mtx

mpif90 -I../../include -std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow main.f90 -L../../lib -lmonolis -lmetis

echo ""
echo "../../bin/gedatsu_simple_mesh_partitioner -n " $1

../../bin/gedatsu_simple_mesh_partitioner -n $1

echo ""
echo "mpirun ./a.out" $1

mpirun -np $1 ./a.out

