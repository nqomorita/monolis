#!/bin/bash

INP=driver/input
OUT=driver/output

./monolis_c_test

mpirun -np 2 ./monolis_c_test

cd input.c

./run.sh
