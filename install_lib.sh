#!/bin/bash

git submodule update --init --recursive
BASE_DIR=$(pwd)

for OPT in "$@"
do
    case $OPT in
        MUMPS)
            FLAG_MUMPS=1
            ;;
        METIS)
            FLAG_METIS=1
            ;;
        SCALAPACK)
            FLAG_SCALAPACK=1
            ;;
    esac
done

if [ "$#" == 0 ]; then
    FLAG_MUMPS=1
    FLAG_METIS=1
    FLAG_SCALAPACK=1
fi

if [ "$FLAG_SCALAPACK" ]; then
    cd submodule/scalapack
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=$BASE_DIR -DCMAKE_C_FLAGS="-Wno-implicit-function-declaration" ..
    make -j
    make install
    cd ../../..
fi

if [ "$FLAG_MUMPS" ]; then
    # under CeCILL license
    cd submodule
    wget http://www.kz.tsukuba.ac.jp/~nmorita/MUMPS_5.3.3.tar.gz
    tar xvf MUMPS_5.3.3.tar.gz
    mv MUMPS_5.3.3 mumps2
    cd mumps2
    #cp Make.inc/Makefile.inc.generic Makefile.inc
    cp ../Makefile.inc.mumps ./
    make d -j
    cp include/*.h ${BASE_DIR}/include
    cp lib/*.a ${BASE_DIR}/lib
#    cd submodule/mumps
#    mkdir build
#    cd build
#    cmake -Darith=d -DCMAKE_INSTALL_PREFIX=$BASE_DIR ..
#    make -j
#    make install
    cd ../../..
fi

if [ "$FLAG_METIS" ]; then
    cd submodule/METIS
    make config prefix=$BASE_DIR
    make install
    cd ../..
fi

