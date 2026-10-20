#!/bin/bash
# monolis GPU ビルド (mog1, NVIDIA HPC SDK 26.1)
# 使用法: ./build_gpu_mog1.sh [deps|lib]
#   deps : 外部ライブラリ (METIS/ParMETIS) とサブモジュールを含めて全ビルド
#   lib  : monolis 本体 + テストのみ再ビルド (既定)
#
# テスト実行時の注意:
#   - 実行バイナリは hpcx OpenMPI の初期化のため mpirun 経由で起動すること
#     (例: mpirun -np 1 ./monolis_test)
#   - MPI 並列 (np >= 2) では GPU 常駐データの袖通信がデバイスポインタで行われるため、
#     CUDA 対応 UCX トランスポートを明示すること:
#       export LD_LIBRARY_PATH=$NVHPC/cuda/lib64:$NVHPC/compilers/lib:$LD_LIBRARY_PATH
#       export UCX_TLS=self,sm,cuda_copy,cuda_ipc
#       mpirun -np 2 -x LD_LIBRARY_PATH -x UCX_TLS ./monolis_test
set -x
export NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.1
export PATH=$NVHPC/comm_libs/mpi/bin:$NVHPC/compilers/bin:$PATH
BASE_DIR=$(cd $(dirname $0) && pwd)
cd $BASE_DIR

MODE=${1:-lib}

if [ "$MODE" = "deps" ]; then
  #> METIS (nvc)
  cd submodule/METIS
  make config prefix=$BASE_DIR cc=nvc
  make install || exit 1
  cd ../..

  #> ParMETIS (mpicc -> nvc)
  cd submodule/ParMETIS
  make config prefix=$BASE_DIR cc=mpicc
  make install || exit 1
  cd ../..

  #> monolis_utils
  cd submodule/monolis_utils
  make clean
  make FLAGS=GPU_GNU || exit 1
  cd ../..

  #> gedatsu
  cd submodule/gedatsu
  make clean
  make FLAGS=SUBMODULE,GPU_GNU || exit 1
  cd ../..

  #> ggtools（driver バイナリは既存ソースの interface 不一致で失敗するが、
  #>          monolis は libggtools.a をリンクしないため失敗を許容する）
  cd submodule/ggtools
  make clean
  make FLAGS=SUBMODULE,GPU_GNU || echo "WARN: ggtools driver build failed (lib build is checked below)"
  test -f lib/libggtools.a || exit 1
  cd ../..
fi

#> monolis (本体 + テストバイナリ)
make all FLAGS=GPU_GNU || exit 1
echo "BUILD_OK"
