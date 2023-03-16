#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_mpi_util_c.h"

void monolis_global_initialize()
{
  monolis_mpi_initialize();
}

void monolis_global_finalize()
{
  monolis_mpi_finalize();
}

void monolis_initialize(
  MONOLIS* mat)
{
  monolis_prm_initialize(mat->PRM);
  monolis_com_initialize(mat->COM);
  monolis_mat_initialize(mat->MAT);
  monolis_mat_initialize(mat->PREC);
}

void monolis_finalize(
  MONOLIS* mat)
{
  monolis_prm_finalize(mat->PRM);
  monolis_com_finalize(mat->COM);
  monolis_mat_finalize(mat->MAT);
  monolis_mat_finalize(mat->PREC);
}
