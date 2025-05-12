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
 * @brief monolis ライブラリの初期化処理
 * @details MPI の初期化関数が呼ばれる。
 * @ingroup def_init
 */
void monolis_global_initialize();

/**
 * @brief monolis ライブラリの終了処理処理
 * @details MPI の終了処理関数が呼ばれる。
 * @ingroup def_init
 */
void monolis_global_finalize();

/**
 * @brief monolis 構造体の初期化処理
 * @param[out] mat monolis 構造体
 * @details 未初期化・終了処理後の monolis 構造体ポインタを入力する。
 * @details 構造体メンバのポインタを NULL 初期化する。構造体メンバのメモリ解放はしない。
 * @ingroup def_init
 */
void monolis_initialize(
  MONOLIS*    mat);

/**
 * @brief monoils 構造体の終了処理
 * @param[inout] mat monolis 構造体
 * @details 初期化処理後の monolis 構造体ポインタを入力する。
 * @details 構造体メンバをメモリ解放し、ポインタ変数を NULL 初期化する。
 * @ingroup def_init
 */
void monolis_finalize(
  MONOLIS* mat);

#ifdef __cplusplus
}
#endif

#endif
