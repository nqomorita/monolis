#!/bin/bash

git submodule update --init --recursive
BASE_DIR=$(pwd)

for OPT in "$@"
do
    case $OPT in
        MUMPS)
            FLAG_MUMPS=1
            FLAG_SCALAPACK=1
            ;;
        SCALAPACK)
            FLAG_SCALAPACK=1
            ;;
        METIS)
            FLAG_METIS=1
            ;;
        METIS64)
            FLAG_METIS64=1
            ;;
    esac
done

if [ "$#" == 0 ]; then
    FLAG_METIS=1
    FLAG_METIS64=0
    FLAG_MUMPS=0
    FLAG_SCALAPACK=0
fi

if [ "$FLAG_SCALAPACK" == 1 ]; then
    cd submodule/scalapack
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=$BASE_DIR -DCMAKE_C_FLAGS="-Wno-implicit-function-declaration" ..
    make -j
    make install
    cd ../../..
fi

if [ "$FLAG_MUMPS" == 1 ]; then
    #cd submodule/mumps
    #mkdir build
    #cd build
    #cmake ..
    #make -DCMAKE_INSTALL_PREFIX=$BASE_DIR -j
    #make install
    #cd ../../..

    # under CeCILL license
    cd submodule
    wget http://www.kz.tsukuba.ac.jp/~nmorita/MUMPS_5.3.3.tar.gz
    tar xvf MUMPS_5.3.3.tar.gz
    mv MUMPS_5.3.3 mumps2
    cd mumps2
    cp ../Makefile.inc.mumps ./Makefile.inc
    make d -j
    cp include/*.h ${BASE_DIR}/include
    cp lib/*.a ${BASE_DIR}/lib
    cd ../..
fi

if [ "$FLAG_METIS" == 1 ]; then
    cd submodule/METIS
    make config prefix=$BASE_DIR
    make install
    cd ../..
fi

if [ "$FLAG_METIS64" == 1 ]; then
    cd submodule/METIS
    make config i64=1 prefix=$BASE_DIR
    make install
    cd ../..
fi
