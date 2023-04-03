#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "monolis_utils.h"
#include "monolis_def_struc_c.h"

void monolis_update_R(
  MONOLIS* mat,
  int      n_dof,
  double*  x)
{
  monolis_mpi_update_R(&mat->com, mat->com.n_internal_vertex, n_dof, x);
}

void monolis_update_I(
  MONOLIS* mat,
  int      n_dof,
  int*     x)
{
  monolis_mpi_update_I(&mat->com, mat->com.n_internal_vertex, n_dof, x);
}

void monolis_update_C(
  MONOLIS*        mat,
  int             n_dof,
  double complex* x)
{
  monolis_mpi_update_C(&mat->com, mat->com.n_internal_vertex, n_dof, x);
}

