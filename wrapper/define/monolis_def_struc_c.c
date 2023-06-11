#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_def_struc_c.h"

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
  monolis_prm_initialize(&mat->prm);
  monolis_mat_initialize(&mat->mat);
  monolis_mat_initialize(&mat->prec);
}

void monolis_finalize(
  MONOLIS* mat)
{
  monolis_prm_initialize(&mat->prm);
  monolis_mat_initialize(&mat->mat);
  monolis_mat_initialize(&mat->prec);
}
