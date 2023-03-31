#!/bin/bash

mpif90 -I../../include -I../../submodule/monolis_utils/include -I../../submodule/gedatsu/include \
-std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o mesher mesher.f90 \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis

./mesher -i mtx.dat

../../submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n 2

mpif90 -I../../include -I../../submodule/monolis_utils/include -I../../submodule/gedatsu/include \
-std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow \
-o solver main.f90 \
-L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis

mpirun -np 2 solver
