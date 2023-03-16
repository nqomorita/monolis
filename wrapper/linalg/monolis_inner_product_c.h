/* monolis_inner_product.h */
#ifndef MONOLIS_INNER_PRODUCT_H
#define MONOLIS_INNER_PRODUCT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include "monolis_def_struc_c.h"

void monolis_inner_product_I(
  MONOLIS* mat,
  int      ndof,
  int*     x,
  int*     y,
  int*     sum);

void monolis_inner_product_R(
  MONOLIS* mat,
  int      ndof,
  double*  x,
  double*  y,
  double*  sum);

void monolis_inner_product_C(
  MONOLIS*        mat,
  int             ndof,
  double complex* x,
  double complex* y,
  double complex* sum);

void monolis_inner_product_I_c_main(
  int  n,
  int  ndof,
  int* x,
  int* y,
  int* sum,
  int  comm);

void monolis_inner_product_R_c_main(
  int     n,
  int     ndof,
  double* x,
  double* y,
  double* sum,
  int     comm);

void monolis_inner_product_C_c_main(
  int             n,
  int             ndof,
  double complex* x,
  double complex* y,
  double complex* sum,
  int             comm);

#ifdef __cplusplus
}
#endif

#endif
