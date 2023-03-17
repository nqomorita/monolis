/* monolis_spmat_handler.h */
#ifndef MONOLIS_SOLVER_H
#define MONOLIS_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>

void monolis_solve_R(
  MONOLIS* mat,
  double*  b,
  double*  x);

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

#ifdef __cplusplus
}
#endif

#endif
