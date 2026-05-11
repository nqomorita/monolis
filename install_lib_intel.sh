#!/bin/bash
set -e

#should be loaded if SQUID
#module load BaseCPU/2023

# 必要な submodule のみを非再帰で init する。
# (--recursive を使うと gedatsu/submodule/monolis_utils などが多重クローン
#  されるため、本プロジェクトで必要なパスを明示的に列挙する)
# Intel 環境では MUMPS/scalapack は別途用意する想定でサブモジュールから除外。
PATHS=(
  submodule/monolis_utils
  submodule/ggtools
  submodule/gedatsu
  submodule/METIS
  submodule/ParMETIS
)

for p in "${PATHS[@]}"; do
  git submodule update --init "$p"
done

# METIS のネスト submodule (GKlib) は METIS リポジトリ内で init する
( cd submodule/METIS && git submodule update --init GKlib )

BASE_DIR=$(pwd)

#> METIS
cd submodule/METIS
make config cc=icx prefix=$BASE_DIR
make install
cd ../..

#> ParMETIS
cd submodule/ParMETIS
make config cc=mpiicx prefix=$BASE_DIR
#make config i64=1 prefix=$BASE_DIR
make install
cd ../..

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
make FLAGS=INTEL
cd ../..

#> gedatsu
cd submodule/gedatsu/
make clean
make FLAGS=INTEL,SUBMODULE
cd ../..

#> ggtools
cd submodule/ggtools/
make clean
make FLAGS=INTEL,SUBMODULE
cd ../..
