#!/bin/bash
set -e

# 3 次元構造格子 SPD 行列を CG + DIAG で MPI 並列 + GPU で解くサンプル
#
# 使い方:
#   ./run.sh                 # デフォルト N=20, NDOF=1, NP=2
#   ./run.sh 40 3 2          # N=40, NDOF=3, NP=2 プロセス

N=${1:-20}
NDOF=${2:-1}
NP=${3:-2}
NLOOP=${4:-1}

# HPC SDK 同梱 MPI (nvfortran ラッパー) を優先
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/26.1/comm_libs/mpi/bin:$PATH
export HWLOC_HIDE_ERRORS=1

# GPU-aware MPI (UCX) 設定:
#   - UCX_MEMTYPE_CACHE=n : デバイスポインタ誤検出によるセグフォを回避
#   - UCX_TLS             : CUDA トランスポートを明示 (sm の host memcpy 誤用を防ぐ)
export UCX_MEMTYPE_CACHE=n
export UCX_TLS=self,sm,cuda_copy,cuda_ipc

# メッシャ (単一プロセス) のビルドと実行 -> node.dat / elem.dat
mpif90 -I../../include -acc -gpu=mem:separate -Mscalapack \
  -o mesher mesher.f90 \
  -L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis

mpirun -np 1 ./mesher -n ${N}

# メッシュを NP 分割 -> parted.0/
mpirun -np 1 ../../bin/gedatsu_simple_mesh_partitioner -n ${NP}

# 並列ソルバ (GPU) のビルドと実行
mpif90 -I../../include -acc -gpu=mem:separate -Mscalapack \
  -o solver main.f90 \
  -L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis

echo "=== single process (np=1) ==="
mpirun -np 1 ./solver -ndof ${NDOF} -nloop ${NLOOP}

echo "=== MPI parallel (np=${NP}) + GPU ==="
mpirun -np ${NP} ./solver -ndof ${NDOF} -nloop ${NLOOP}
