/* monolis_spmat_handler.h */
#ifndef MONOLIS_SOLVE_H
#define MONOLIS_SOLVE_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 線形ソルバ関数（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] com 通信テーブル構造体
 * @param[in] b 右辺ベクトル
 * @param[inout] x 解ベクトル
 * @ingroup solver
 */
void monolis_solve_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  double*      b,
  double*      x);

/**
 * @brief 線形ソルバ関数（メイン関数）
 * @param[in] n
 * @param[in] np
 * @param[in] nz
 * @param[in] n_dof
 * @param[in] A
 * @param[in] x
 * @param[in] b
 * @param[in] index
 * @param[in] item
 * @param[in] my_rank
 * @param[in] comm
 * @param[in] comm_size
 * @param[in] recv_n_neib
 * @param[in] recv_nitem
 * @param[in] recv_neib_pe
 * @param[in] recv_index
 * @param[in] recv_item
 * @param[in] send_n_neib
 * @param[in] send_nitem
 * @param[in] send_neib_pe
 * @param[in] send_index
 * @param[in] send_item
 * @param[in] Iarray
 * @param[in] Rarray
 * @ingroup dev_solver
 */
void monolis_solve_R_c_main(
  int     n,
  int     np,
  int     nz,
  int     n_dof,
  double* A,
  double* x,
  double* b,
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
  int*    send_item,
  int*    Iarray,
  double* Rarray);

/**
 * @brief 線形ソルバ関数（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] com 通信テーブル構造体
 * @param[in] b 右辺ベクトル
 * @param[inout] x 解ベクトル
 * @ingroup solver
 */
void monolis_solve_C(
  MONOLIS*        mat,
  MONOLIS_COM*    com,
  double _Complex* b,
  double _Complex* x);

/**
 * @brief 線形ソルバ関数（メイン関数）
 * @param[in] n
 * @param[in] np
 * @param[in] nz
 * @param[in] n_dof
 * @param[in] A
 * @param[in] x
 * @param[in] b
 * @param[in] index
 * @param[in] item
 * @param[in] my_rank
 * @param[in] comm
 * @param[in] comm_size
 * @param[in] recv_n_neib
 * @param[in] recv_nitem
 * @param[in] recv_neib_pe
 * @param[in] recv_index
 * @param[in] recv_item
 * @param[in] send_n_neib
 * @param[in] send_nitem
 * @param[in] send_neib_pe
 * @param[in] send_index
 * @param[in] send_item
 * @param[in] Iarray
 * @param[in] Rarray
 * @ingroup dev_solver
 */
void monolis_solve_C_c_main(
  int     n,
  int     np,
  int     nz,
  int     n_dof,
  double _Complex* A,
  double _Complex* x,
  double _Complex* b,
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
  int*    send_item,
  int*    Iarray,
  double* Rarray);

#ifdef __cplusplus
}
#endif

#endif
