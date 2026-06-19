#!/bin/bash
set -e

# CVAE サンプル: クラスラベルを条件として 1D 波形データを学習し、
# 条件付き再構成・条件付き生成を行う
#
# 使い方:
#   ./run.sh

mpif90 -I../../include \
  -std=legacy \
  -o solver main.f90 \
  -L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils \
  -lmetis -llapack -lblas

./solver
