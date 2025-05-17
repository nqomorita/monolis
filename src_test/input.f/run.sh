#!/bin/bash

mpif90 -I../../include \
-std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o mesher mesher.f90 \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -llapack -lblas

mpirun -np 1 ./mesher -i mtx.dat

mpirun -np 1 ../../bin/gedatsu_simple_mesh_partitioner -n 3

mpif90 -I../../include \
-std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o solver main.f90 \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -lscalapack -llapack -lblas

#mpirun -np 1 ./solver

#mpirun -np 3 ./solver

mpif90 -I../../include \
-std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o solver_arbit main.arbit.f90 \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -lscalapack -llapack -lblas

mpirun -np 1 ./solver_arbit

#mpirun -np 3 ./solver_arbit
