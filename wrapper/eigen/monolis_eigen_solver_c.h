/* monolis_eigen_solver.h */
#ifndef MONOLIS_EIGEN_SOLVER_H
#define MONOLIS_EIGEN_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Lanczos 法（逆反復、最小固有値、実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] com 通信テーブル構造体
 * @param[inout] n_get_eigen 取得固有値数
 * @param[in] ths 収束判定閾値
 * @param[in] maxiter 最大反復回数
 * @param[inout] eigen_value 固有値
 * @param[inout] eigen_mode 固有ベクトル
 * @param[in] is_Dirichlet_bc Dirhchlet 境界条件判定フラグ
 * @ingroup eigen
 */
void monolis_eigen_inverted_standard_lanczos_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  int*         n_get_eigen,
  double       ths,
  int          maxiter,
  double*      eigen_value,
  double**     eigen_mode,
  bool*        is_Dirichlet_bc);

/**
 * @brief !
 * @param[in] n
 * @param[in] np
 * @param[in] nz
 * @param[in] n_dof
 * @param[in] A
 * @param[in] index
 * @param[in] item
 * @param[in] my_rank
 * @param[in] comm
 * @param[in] comm_size
 * @param[in] recv_n_neib
 * @param[in] recv_nitem
 * @param[in] revc_neib_pe
 * @param[in] recv_index
 * @param[in] recv_item
 * @param[in] send_n_neib
 * @param[in] send_nitem
 * @param[in] send_neib_pe
 * @param[in] send_index
 * @param[in] send_item
 * @param[in] Iarray
 * @param[in] Rarray
 * @param[in] n_get_eigen
 * @param[in] ths
 * @param[in] eigen_maxiter
 * @param[in] eigen_value
 * @param[in] eigen_mode_tmp
 * @param[in] is_Dirichlet_bc_int
 * @ingroup dev_solver
 */
void monolis_eigen_inverted_standard_lanczos_R_c_main(
  int     n,
  int     np,
  int     nz,
  int     n_dof,
  double* A,
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
  double* Rarray,
  int*    n_get_eigen,
  double  ths,
  int     eigen_maxiter,
  double* eigen_value,
  double* eigen_mode_tmp,
  int*    is_Dirichlet_bc_int);

/**
 * @brief Lanczos 法（順反復、最大固有値、実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] com 通信テーブル構造体
 * @param[inout] n_get_eigen 取得固有値数
 * @param[in] ths 収束判定閾値
 * @param[in] maxiter 最大反復回数
 * @param[inout] eigen_value 固有値
 * @param[inout] eigen_mode 固有ベクトル
 * @param[in] is_Dirichlet_bc Dirhchlet 境界条件判定フラグ
 * @ingroup eigen
 */
void monolis_eigen_standard_lanczos_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  int*         n_get_eigen,
  double       ths,
  int          maxiter,
  double*      eigen_value,
  double**     eigen_mode,
  bool*        is_Dirichlet_bc);

/**
 * @brief !
 * @param[in] n
 * @param[in] np
 * @param[in] nz
 * @param[in] n_dof
 * @param[in] A
 * @param[in] index
 * @param[in] item
 * @param[in] my_rank
 * @param[in] comm
 * @param[in] comm_size
 * @param[in] recv_n_neib
 * @param[in] recv_nitem
 * @param[in] revc_neib_pe
 * @param[in] recv_index
 * @param[in] recv_item
 * @param[in] send_n_neib
 * @param[in] send_nitem
 * @param[in] send_neib_pe
 * @param[in] send_index
 * @param[in] send_item
 * @param[in] Iarray
 * @param[in] Rarray
 * @param[in] n_get_eigen
 * @param[in] ths
 * @param[in] eigen_maxiter
 * @param[in] eigen_value
 * @param[in] eigen_mode_tmp
 * @param[in] is_Dirichlet_bc_int
 * @ingroup dev_solver
 */
void monolis_eigen_standard_lanczos_R_c_main(
  int     n,
  int     np,
  int     nz,
  int     n_dof,
  double* A,
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
  double* Rarray,
  int*    n_get_eigen,
  double  ths,
  int     eigen_maxiter,
  double* eigen_value,
  double* eigen_mode_tmp,
  int*    is_Dirichlet_bc_int);

#ifdef __cplusplus
}
#endif

#endif
