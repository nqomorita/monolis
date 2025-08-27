#!/bin/bash

#should be loaded
#module load fj
#module load fjmpi
#module load metis
#module load parmetis

git submodule update --init --recursive
BASE_DIR=$(pwd)

#> METIS
#cd submodule/METIS
#make config cc=fccpx prefix=$BASE_DIR
#make install
#cd ../..

#> MUMPS
#cd submodule/mumps
#mkdir build
#cd build
#cmake ..
#make -DCMAKE_INSTALL_PREFIX=$BASE_DIR -j
#cd ../../..

#> monolis_utils
cd submodule/monolis_utils/
make clean
make FLAGS=A64FX
cd ../..

#> gedatsu
cd submodule/gedatsu/
make clean
make FLAGS=A64FX,SUBMODULE
cd ../..

#> gedatsu
cd submodule/ggtools/
make clean
make FLAGS=A64FX,SUBMODULE
cd ../..
