#!/bin/bash

make clean
make

./mesher.out -i ./bcsstk16.mtx

#../../bin/gedatsu_simple_mesh_partitioner -n $1

./a.out
