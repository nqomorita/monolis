#!/bin/bash

mpif90 -I../../include -I../../submodule/monolis_utils/include -std=legacy -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow mesher.f90 \
-L../../lib -lmonolis -L../../submodule/monolis_utils/lib -lmonolis_utils -lmetis

./a.out -i mtx.dat

../../submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n 2 -o elem.dat
