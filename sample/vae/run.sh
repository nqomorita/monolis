#!/bin/bash
set -e

# VAE サンプル: 合成 1D 波形データを学習し、再構成・サンプリングを行う
#
# 使い方:
#   ./run.sh

mpif90 -I../../include \
  -std=legacy \
  -o solver main.f90 \
  -L../../lib -lmonolis_solver -lgedatsu -lmonolis_utils \
  -lmetis -llapack -lblas

./solver
