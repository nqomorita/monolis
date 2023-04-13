/* monolis_def_mat.h */
#ifndef MONOLIS_MATVEC_H
#define MONOLIS_MATVEC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include "monolis_def_struc_c.h"

void monolis_matvec_product_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  double*      x,
  double*      y);

void monolis_matvec_product_C(
  MONOLIS*        mat,
  MONOLIS_COM*    com,
  double complex* x,
  double complex* y);

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
  double complex* A,
  double complex* x,
  double complex* y,
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
