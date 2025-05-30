#!/bin/bash

echo "mesher"

mpif90 -I../../include \
-std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o mesher mesher.f90 \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -llapack -lblas

echo "solver C"

mpicc -I../../include \
-c -o main.o main.c

mpif90 -I../../include -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o csolver main.o \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -llapack -lblas

for A in 100 1000 10000
do
  for B in 1 0.01 0.0001 0.000001 0.00000001 0.0000000001 0.000000000001
  do
    for C in 1 2 4 8 16
    do
      echo " ###ILL" $A $B $C
      mpirun -np 1 ./mesher -n $A -g $B
      mpirun -np 1 ../../bin/gedatsu_simple_mesh_partitioner -in node.c.dat -ie elem.c.dat -n $C
      mpirun --oversubscribe -np $C ./csolver
    done
  done
done

