#!/bin/bash

./monolis_test

mpirun -np 2 ./monolis_test

cd input.f

./run.sh
