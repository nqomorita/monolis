/* monolis_spmat_handler.h */
#ifndef MONOLIS_SOLVE_H
#define MONOLIS_SOLVE_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup solver
 */
void monolis_solve_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  double*      b,
  double*      x);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
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
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup solver
 */
void monolis_solve_C(
  MONOLIS*        mat,
  MONOLIS_COM*    com,
  double _Complex* b,
  double _Complex* x);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
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
