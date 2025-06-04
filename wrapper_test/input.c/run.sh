#!/bin/bash

echo "mesher"

mpif90 -I../../include \
-std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o mesher mesher.f90 \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -llapack -lblas

mpirun -np 1 ./mesher -i mtx.dat

../../bin/gedatsu_simple_mesh_partitioner -n 3

echo "solver"

mpicc -I../../include \
-c -o main.o main.c

mpif90 -I../../include \
-o solver main.o \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -lscalapack -llapack -lblas

echo "serial"

mpirun -np 1 ./solver

echo "parallel"

mpirun -np 3 ./solver

echo "solver arbit"

mpicc -I../../include \
-c -o main.arbit.o main.arbit.c

mpif90 -I../../include \
-o solver.arbit main.arbit.o \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -lscalapack -llapack -lblas

echo "serial"

mpirun -np 1 ./solver.arbit

echo "parallel"

mpirun -np 3 solver.arbit

