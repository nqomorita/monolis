#!/bin/bash
set -e

# 必要な submodule のみを非再帰で init する。
# (--recursive を使うと gedatsu/submodule/monolis_utils などが多重クローン
#  されるため、本プロジェクトで必要なパスを明示的に列挙する)
PATHS=(
  submodule/monolis_utils
  submodule/ggtools
  submodule/gedatsu
  submodule/METIS
  submodule/ParMETIS
  submodule/mumps
  submodule/scalapack
)

for p in "${PATHS[@]}"; do
  git submodule update --init "$p"
done

# METIS のネスト submodule (GKlib) は METIS リポジトリ内で init する
( cd submodule/METIS && git submodule update --init GKlib )

BASE_DIR=$(pwd)

#> scalapack
#cd submodule/scalapack
#mkdir build
#cd build
#export CMAKE_POLICY_VERSION_MINIMUM=3.5
#cmake -DCMAKE_INSTALL_PREFIX=$BASE_DIR -DCMAKE_C_FLAGS="-Wno-implicit-function-declaration" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DSCALAPACK_BUILD_TESTS=OFF ..
#make -j
#make install
#unset CMAKE_POLICY_VERSION_MINIMUM
#cd ../../..

#> METIS
cd submodule/METIS
make config prefix=$BASE_DIR cc=nvc
make install
cd ../..
#make config i64=1 prefix=$BASE_DIR

#> ParMETIS
cd submodule/ParMETIS
make config prefix=$BASE_DIR
#make config i64=1 prefix=$BASE_DIR
make install
cd ../..

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
make FLAGS=GPU_GNU lib
cd ../..

#> gedatsu
cd submodule/gedatsu/
make clean
make FLAGS=SUBMODULE,GPU_GNU
cd ../..

#> ggtools
cd submodule/ggtools/
make clean
make FLAGS=SUBMODULE,GPU_GNU
cd ../..
