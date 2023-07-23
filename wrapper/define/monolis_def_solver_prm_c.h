/* monolis_def_struc.h */
#ifndef MONOLIS_DEF_SOLVER_H
#define MONOLIS_DEF_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

#include "monolis_utils.h"

static const int MONOLIS_PRM_IARRAY_SIZE = 100;

static const int MONOLIS_PRM_RARRAY_SIZE = 100;

static const int MONOLIS_ITER_CG       = 1;

static const int MONOLIS_ITER_GROPCG   = 2;

static const int MONOLIS_ITER_PIPECG   = 3;

static const int MONOLIS_ITER_PIPECR   = 4;

static const int MONOLIS_ITER_BICGSTAB = 5;

static const int MONOLIS_ITER_PIPEBICGSTAB = 6;

static const int MONOLIS_ITER_BICGSTAB_NOPREC = 7;

//static const int MONOLIS_ITER_CABICGSTAB_NOPREC = 8;

static const int MONOLIS_ITER_PIPEBICGSTAB_NOPREC = 8;

//static const int MONOLIS_ITER_GMRES = 10;

static const int MONOLIS_ITER_COCG = 9;

static const int MONOLIS_PREC_NONE   = 0;

static const int MONOLIS_PREC_DIAG   = 1;

//static const int MONOLIS_PREC_ILU    = 2;

//static const int MONOLIS_PREC_JACOBI = 3;

static const int MONOLIS_PREC_SOR    = 2;

//static const int MONOLIS_PREC_SPIKE  = 5;

//static const int MONOLIS_PREC_LU     = 6;

//static const int MONOLIS_PREC_MUMPS  = 7;

//static const int MONOLIS_PREC_ROM    = 8;

//static const int MONOLIS_PREC_MF     = 9;

//static const int MONOLIS_PREC_MUMPS_LOCAL = 10;

static const int MONOLIS_PRM_I_METHOD = 1;

static const int MONOLIS_PRM_I_PRECOND = 2;

static const int MONOLIS_PRM_I_MAX_ITER = 3;

static const int MONOLIS_PRM_I_CUR_ITER = 4;

static const int MONOLIS_PRM_I_IERR = 5;

//static const int MONOLIS_PRM_I_IS_SCALING = 6;

//static const int MONOLIS_PRM_I_IS_REORDERING = 7;

static const int MONOLIS_PRM_I_IS_INIT_X = 8;

static const int MONOLIS_PRM_I_IS_SYM_MATRIX = 9;

static const int MONOLIS_PRM_I_IS_DEBUG = 10;

static const int MONOLIS_PRM_I_IS_MEASUREMENT = 11;

static const int MONOLIS_PRM_I_IS_CHECK_DIAG = 12;

static const int MONOLIS_PRM_I_IS_PREC_STORED = 13;

static const int MONOLIS_PRM_I_SHOW_ITERLOG = 14;

static const int MONOLIS_PRM_I_SHOW_TIME = 15;

static const int MONOLIS_PRM_I_SHOW_SUMMARY = 16;

static const int MONOLIS_PRM_I_SHOW_TIME_STATISTICS = 17;

static const int MONOLIS_PRM_R_TOL = 1;

static const int MONOLIS_PRM_R_CUR_RESID = 2;

static const int MONOLIS_R_TIME_SOL = 3;

static const int MONOLIS_R_TIME_PREP = 4;

static const int MONOLIS_R_TIME_SPMV = 5;

static const int MONOLIS_R_TIME_DOTP = 6;

static const int MONOLIS_R_TIME_PREC = 7;

static const int MONOLIS_R_TIME_COMM_DOTP = 8;

static const int MONOLIS_R_TIME_COMM_SPMV = 9;

typedef struct {
  char   com_top_dir_name[1024];
  char   com_part_dir_name[1024];
  char   com_file_name[1024];
  int    Iarray[100];
  double Rarray[100];
} MONOLIS_PRM;

/**
 * @brief パラメータ 構造体の初期化処理
 * @param[out] prm パラメータ構造体
 * @ingroup def_init
 */
void monolis_prm_initialize(
  MONOLIS_PRM* prm);

/**
 * @brief パラメータ構造体の終了処理
 * @param[inout] prm パラメータ構造体
 * @ingroup def_init
 */
void monolis_prm_finalize(
  MONOLIS_PRM* prm);

#ifdef __cplusplus
}
#endif

#endif
