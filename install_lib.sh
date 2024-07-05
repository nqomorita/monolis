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
cmake ..
make -DCMAKE_INSTALL_PREFIX=$BASE_DIR -j 
cd ../../..

#> monolis_utils
cd submodule/monolis_utils/
make
cd ../..

#> gedatsu
cd submodule/gedatsu/
make FLAGS=SUBMODULE
cd ../..
