#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "monolis_utils.h"
#include "monolis_inner_product_c.h"
#include "monolis_def_struc_c.h"

void monolis_inner_product_I(
  MONOLIS* mat,
  int      ndof,
  int*     x,
  int*     y,
  int*     sum)
{
  int n = mat->mat.N;
  if(mat->com.comm_size > 1) n = mat->com.n_internal_vertex;
  int comm = mat->com.comm;

  monolis_inner_product_I_c_main(n, ndof, x, y, sum, comm);
}

void monolis_inner_product_R(
  MONOLIS* mat,
  int      ndof,
  double*  x,
  double*  y,
  double*  sum)
{
  int n = mat->mat.N;
  if(mat->com.comm_size > 1) n = mat->com.n_internal_vertex;
  int comm = mat->com.comm;

  monolis_inner_product_R_c_main(n, ndof, x, y, sum, comm);
}

void monolis_inner_product_C(
  MONOLIS*        mat,
  int             ndof,
  double complex* x,
  double complex* y,
  double complex* sum)
{
  int n = mat->mat.N;
  if(mat->com.comm_size > 1) n = mat->com.n_internal_vertex;
  int comm = mat->com.comm;

  monolis_inner_product_C_c_main(n, ndof, x, y, sum, comm);
}
