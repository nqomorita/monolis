#!/bin/bash
set -e

#should be loaded
#module load fj
#module load fjmpi
#module load metis
#module load parmetis

# 必要な submodule のみを非再帰で init する。
# (--recursive を使うと gedatsu/submodule/monolis_utils などが多重クローン
#  されるため、本プロジェクトで必要なパスを明示的に列挙する)
# A64FX 環境では METIS/ParMETIS/MUMPS/scalapack は module load を想定し、
# サブモジュールとしては取得しない。
PATHS=(
  submodule/monolis_utils
  submodule/ggtools
  submodule/gedatsu
)

for p in "${PATHS[@]}"; do
  git submodule update --init "$p"
done

BASE_DIR=$(pwd)

#> METIS
#cd submodule/METIS
#make config cc=fccpx prefix=$BASE_DIR
#make install
#cd ../..

#> ParMETIS
#cd submodule/ParMETIS
#make config prefix=$BASE_DIR
#make config i64=1 prefix=$BASE_DIR
#make install
#cd ../..

#> MUMPS
#cd submodule/mumps
#mkdir build
#cd build
#cmake -DCMAKE_INSTALL_PREFIX=$BASE_DIR ..
#make -j
#make install
#cd ../../..

#> monolis_utils
cd submodule/monolis_utils/
make clean
make FLAGS=A64FX lib
cd ../..

#> gedatsu
cd submodule/gedatsu/
make clean
make FLAGS=A64FX,SUBMODULE
cd ../..

#> ggtools
cd submodule/ggtools/
make clean
make FLAGS=A64FX,SUBMODULE
cd ../..
