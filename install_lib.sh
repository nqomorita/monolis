#!/bin/bash

git submodule update --init --recursive
BASE_DIR=$(pwd)

#> for METIS
cd submodule/METIS
make config prefix=$BASE_DIR
make install
cd ../..

#> for MUMPS
cd submodule/mumps
mkdir build
cd build
cmake -Darith=d -DCMAKE_INSTALL_PREFIX=$BASE_DIR ..
make -j
make install

