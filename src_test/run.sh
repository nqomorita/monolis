#!/bin/bash

./monolis_test | tee test_list.dat

mpirun -np 2 ./monolis_test | tee -a test_list.dat

cd input.f

./run.sh
