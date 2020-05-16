#!/bin/bash

#> for METIS
BASE_DIR=$(pwd)
git submodule update --init --recursive
cd submodule/METIS
make config prefix=$BASE_DIR
make install
