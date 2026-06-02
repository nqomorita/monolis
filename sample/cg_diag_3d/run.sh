#!/bin/bash
set -e

# 3 次元構造格子 SPD 行列を CG + DIAG で解くサンプル
#
# 使い方:
#   ./run.sh           # デフォルト N=10 NDOF=1 で 1 / 2 プロセス実行
#   ./run.sh 20 3      # N=20, NDOF=3 で実行

N=${1:-10}
NDOF=${2:-1}
NLOOP=${3:-1}

mpif90 -I../../include -acc \
  -o solver main.f90 \
  -L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils \
  -lmetis 

mpirun -np 1 ./solver -n ${N} -ndof ${NDOF} -nloop ${NLOOP}


