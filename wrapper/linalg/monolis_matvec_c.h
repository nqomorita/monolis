/* monolis_def_mat.h */
#ifndef MONOLIS_MATVEC_H
#define MONOLIS_MATVEC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "monolis_def_struc_c.h"

/**
 * @brief 疎行列ベクトル積（実数型）
 * @param[in] mat monolis 構造体
 * @param[in] com COM 構造体
 * @param[inout] x 右辺ベクトル
 * @param[out] y 結果ベクトル
 * @ingroup linalg
 */
void monolis_matvec_product_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  double*      x,
  double*      y);

/**
 * @brief 疎行列ベクトル積（複素数型）
 * @param[in] mat monolis 構造体
 * @param[in] com COM 構造体
 * @param[inout] x 右辺ベクトル
 * @param[out] y 結果ベクトル
 * @ingroup linalg
 */
void monolis_matvec_product_C(
  MONOLIS*        mat,
  MONOLIS_COM*    com,
  double _Complex* x,
  double _Complex* y);

void monolis_matvec_product_R_c_main(
  int     n,
  int     np,
  int     nz,
  int     ndof,
  double* A,
  double* x,
  double* y,
  int*    index,
  int*    item,
  int     my_rank,
  int     comm,
  int     comm_size,
  int     recv_n_neib,
  int     recv_nitem,
  int*    recv_neib_pe,
  int*    recv_index,
  int*    recv_item,
  int     send_n_neib,
  int     send_nitem,
  int*    send_neib_pe,
  int*    send_index,
  int*    send_item);

void monolis_matvec_product_C_c_main(
  int             n,
  int             np,
  int             nz,
  int             ndof,
  double _Complex* A,
  double _Complex* x,
  double _Complex* y,
  int*            index,
  int*            item,
  int             my_rank,
  int             comm,
  int             comm_size,
  int             recv_n_neib,
  int             recv_nitem,
  int*            recv_neib_pe,
  int*            recv_index,
  int*            recv_item,
  int             send_n_neib,
  int             send_nitem,
  int*            send_neib_pe,
  int*            send_index,
  int*            send_item);

#ifdef __cplusplus
}
#endif

#endif
