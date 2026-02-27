#!/bin/bash

make clean
make

#./mesher.out -i ./input/A.mtx
#./mesher.out -i ./test.mtx
#./mesher.out -i ./test.3.mtx
./mesher.out -i ./bcsstk16.mtx
#./mesher.out -i ./s3dkq4m2.mtx
#./mesher.out -i ./input/beam.mtx

#../../bin/gedatsu_simple_mesh_partitioner -n $1

./a.out

