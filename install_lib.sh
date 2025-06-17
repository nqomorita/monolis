#!/bin/bash

git submodule update --init --recursive
BASE_DIR=$(pwd)

#> scalapack
cd submodule/scalapack
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$BASE_DIR -DCMAKE_C_FLAGS="-Wno-implicit-function-declaration" ..
make -j
make install
cd ../../..

#> METIS
cd submodule/METIS
make config prefix=$BASE_DIR
make install
cd ../..
#make config i64=1 prefix=$BASE_DIR

#> MUMPS
cd submodule/mumps
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$BASE_DIR ..
make -j 
make install
cd ../../..

#> monolis_utils
cd submodule/monolis_utils/
make clean
make
cd ../..

#> gedatsu
cd submodule/gedatsu/
make clean
make FLAGS=SUBMODULE
cd ../..

#> gedatsu
cd submodule/ggtools/
make clean
make FLAGS=SUBMODULE
cd ../..
