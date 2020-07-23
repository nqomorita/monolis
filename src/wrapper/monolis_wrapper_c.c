#include <stdio.h>
#include <stdlib.h>
#include "monolis.h"

void monolis_initialize(
  MONOLIS *mat)
{
  printf("monolis_initialize in C\n");
}

void monolis_finalize(
  MONOLIS *mat)
{
  printf("monolis_finalize in C\n");
}

void monolis_get_nonzero_pattern(
  MONOLIS* mat,
  int      nnode,
  int      nbase_func,
  int      ndof,
  int      nelem,
  int**    elem)
{
  printf("monolis_get_nonzero_pattern in C\n");
}

void monolis_add_sparse_matrix(
  MONOLIS* mat,
  int      nbase_func,
  int      *connectivity,
  double** local_mat)
{
  printf("monolis_add_sparse_matrix in C\n");
}

void monolis_set_Dirichlet_bc(
  MONOLIS* mat,
  double*  b,
  int      node_id,
  int      ndof_bc,
  double*  val)
{
  printf("monolis_set_Dirichlet_bc in C\n");
}

void monolis_solve(
  MONOLIS* mat,
  double*  b,
  double*  x)
{
  printf("monolis_solve in C\n");
}
