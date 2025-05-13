#!/bin/bash

#should be loaded if SQUID
#module load BaseCPU/2023

git submodule update --init --recursive
BASE_DIR=$(pwd)

#> METIS
cd submodule/METIS
make config cc=icx prefix=$BASE_DIR
make install
cd ../..

#> MUMPS
#cd submodule/mumps
#mkdir build
#cd build
#cmake ..
#make -DCMAKE_INSTALL_PREFIX=$BASE_DIR -j
#cd ../../..

#> monolis_utils
cd submodule/monolis_utils/
make FLAGS=INTEL
cd ../..

#> gedatsu
cd submodule/gedatsu/
make FLAGS=INTEL,SUBMODULE
cd ../..
