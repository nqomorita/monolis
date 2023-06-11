/* monolis_inner_product.h */
#ifndef MONOLIS_INNER_PRODUCT_H
#define MONOLIS_INNER_PRODUCT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "monolis_def_struc_c.h"

void monolis_inner_product_I(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  int          n_dof,
  int*         x,
  int*         y,
  int*         sum);

void monolis_inner_product_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  int          n_dof,
  double*      x,
  double*      y,
  double*      sum);

void monolis_inner_product_C(
  MONOLIS*        mat,
  MONOLIS_COM*    com,
  int             n_dof,
  double _Complex* x,
  double _Complex* y,
  double _Complex* sum);

void monolis_inner_product_I_c_main(
  int  n,
  int  n_dof,
  int* x,
  int* y,
  int* sum,
  int  comm);

void monolis_inner_product_R_c_main(
  int     n,
  int     n_dof,
  double* x,
  double* y,
  double* sum,
  int     comm);

void monolis_inner_product_C_c_main(
  int             n,
  int             n_dof,
  double _Complex* x,
  double _Complex* y,
  double _Complex* sum,
  int             comm);

#ifdef __cplusplus
}
#endif

#endif
