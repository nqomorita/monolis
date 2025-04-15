#!/bin/bash

INP=driver/input
OUT=driver/output

mpirun -np 1 ./monolis_c_test

mpirun -np 2 ./monolis_c_test

cd input.c

./run.sh
