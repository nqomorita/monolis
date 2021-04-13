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
    cmake -DCMAKE_INSTALL_PREFIX=$BASE_DIR ..
    make -j
    make install
    cd ../../..
fi

if [ "$FLAG_MUMPS" ]; then
    cd submodule/mumps
    mkdir build
    cd build
    cmake -Darith=d -DCMAKE_INSTALL_PREFIX=$BASE_DIR ..
    make -j
    make install
    cd ../../..
fi

if [ "$FLAG_METIS" ]; then
    cd submodule/METIS
    make config prefix=$BASE_DIR
    make install
    cd ../..
fi
