/* monolis_def_struc.h */
#ifndef MONOLIS_DEF_STRUC_H
#define MONOLIS_DEF_STRUC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "monolis_utils.h"
#include "monolis_def_solver_prm_c.h"
#include "monolis_def_mat_c.h"

typedef struct {
  MONOLIS_PRM prm;
  MONOLIS_MAT mat;
  MONOLIS_MAT prec;
} MONOLIS;

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_init
 */
void monolis_global_initialize();

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_init
 */
void monolis_global_finalize();

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_init
 */
void monolis_initialize(
  MONOLIS*    mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_init
 */
void monolis_finalize(
  MONOLIS* mat);

#ifdef __cplusplus
}
#endif

#endif
