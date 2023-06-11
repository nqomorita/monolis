#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_inner_product_c.h"
#include "monolis_def_struc_c.h"

void monolis_inner_product_I(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  int          ndof,
  int*         x,
  int*         y,
  int*         sum)
{
  int n = mat->mat.N;
  if(com->comm_size > 1) n = com->n_internal_vertex;

  monolis_inner_product_I_c_main(n, ndof, x, y, sum, com->comm);
}

void monolis_inner_product_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  int          ndof,
  double*      x,
  double*      y,
  double*      sum)
{
  int n = mat->mat.N;
  if(com->comm_size > 1) n = com->n_internal_vertex;

  monolis_inner_product_R_c_main(n, ndof, x, y, sum, com->comm);
}

void monolis_inner_product_C(
  MONOLIS*        mat,
  MONOLIS_COM*    com,
  int             ndof,
  double _Complex* x,
  double _Complex* y,
  double _Complex* sum)
{
  int n = mat->mat.N;
  if(com->comm_size > 1) n = com->n_internal_vertex;

  monolis_inner_product_C_c_main(n, ndof, x, y, sum, com->comm);
}
