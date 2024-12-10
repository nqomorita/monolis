#!/bin/bash

echo "mesher"

mpif90 -I../../include \
-std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o mesher mesher.f90 \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -llapack -lblas

mpirun -np 1 ./mesher -n 100 -g 0.1

mpirun -np 1 ../../bin/gedatsu_simple_mesh_partitioner -in node.f.dat -ie elem.f.dat -n 3
mpirun -np 1 ../../bin/gedatsu_simple_mesh_partitioner -in node.c.dat -ie elem.c.dat -n 3


echo "solver F"

mpif90 -I../../include \
-std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o fsolver main.f90 \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -llapack -lblas

mpirun -np 1 ./fsolver

mpirun -np 3 ./fsolver


echo "solver C"

mpicc -I../../include \
-c -o main.o main.c

mpif90 -I../../include -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o csolver main.o \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -llapack -lblas

mpirun -np 1 ./csolver

mpirun -np 3 ./csolver
