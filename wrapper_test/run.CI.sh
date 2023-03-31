#!/bin/bash

INP=driver/input
OUT=driver/output

./monolis_c_test

mpirun --oversubscribe --allow-run-as-root -np 2 ./monolis_c_test

cd input.c

./run.CI.sh
