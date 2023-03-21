#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_def_struc_c.h"

void monolis_com_input_comm_table(
  MONOLIS* mat,
  const char* input_file_dir)
{
  if(mat->com.comm_size <= 1){
    mat->com.comm_size = 1;
    mat->com.send_n_neib = 0;
    mat->com.recv_n_neib = 0;
    mat->com.send_neib_pe = (int*)calloc(1, sizeof(int));
    mat->com.send_index = (int*)calloc(2, sizeof(int));
    mat->com.send_item = (int*)calloc(1, sizeof(int));
    mat->com.recv_neib_pe = (int*)calloc(1, sizeof(int));
    mat->com.recv_index = (int*)calloc(2, sizeof(int));
    mat->com.recv_item = (int*)calloc(1, sizeof(int));
    return;
  }
}

void monolis_global_initialize()
{
  monolis_mpi_initialize();
}

void monolis_global_finalize()
{
  monolis_mpi_finalize();
}

void monolis_initialize(
  MONOLIS* mat,
  const char* input_file_dir)
{
  monolis_prm_initialize(&mat->prm);
  monolis_com_initialize(&mat->com);
  monolis_mat_initialize(&mat->mat);
  monolis_mat_initialize(&mat->prec);
  monolis_com_input_comm_table(mat,input_file_dir);
}

void monolis_finalize(
  MONOLIS* mat)
{
  monolis_prm_initialize(&mat->prm);
  monolis_com_initialize(&mat->com);
  monolis_mat_initialize(&mat->mat);
  monolis_mat_initialize(&mat->prec);
}
