#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_def_struc_c.h"

void monolis_com_input_comm_table(
  MONOLIS* mat,
  const char* top_dir_name,
  const char* part_dir_name,
  const char* file_name)
{
  if(mat->com.comm_size <= 1){
    mat->com.comm_size = 1;
    mat->com.send_n_neib = 0;
    mat->com.recv_n_neib = 0;
    mat->com.send_neib_pe = monolis_alloc_I_1d(mat->com.send_neib_pe, 1);
    mat->com.send_index = monolis_alloc_I_1d(mat->com.send_index, 2);
    mat->com.send_item = monolis_alloc_I_1d(mat->com.send_item, 1);
    mat->com.recv_neib_pe = monolis_alloc_I_1d(mat->com.recv_neib_pe, 1);
    mat->com.recv_index = monolis_alloc_I_1d(mat->com.recv_index, 2);
    mat->com.recv_item = monolis_alloc_I_1d(mat->com.recv_item, 1);
    return;
  }

  const char* input_file_name;
  char* header;
  char* body;

  header = (char*)malloc(sizeof(char)*MONOLIS_CHARLEN);
  snprintf(header, MONOLIS_CHARLEN, "%s/%s", top_dir_name, part_dir_name);

  body = (char*)malloc(sizeof(char)*MONOLIS_CHARLEN);
  snprintf(body, MONOLIS_CHARLEN, "%s.recv", file_name);
  input_file_name = monolis_get_output_file_name_by_domain_id(header, body, mat->com.my_rank);
  monolis_input_recv_com_table(input_file_name, &mat->com);

  snprintf(body, MONOLIS_CHARLEN, "%s.send", file_name);
  input_file_name = monolis_get_output_file_name_by_domain_id(header, body, mat->com.my_rank);
  monolis_input_send_com_table(input_file_name, &mat->com);

  snprintf(body, MONOLIS_CHARLEN, "%s.n_internal", file_name);
  input_file_name = monolis_get_output_file_name_by_domain_id(header, body, mat->com.my_rank);
  monolis_input_internal_vertex_number(input_file_name, &mat->com.n_internal_vertex);
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
  MONOLIS* mat)
{
  monolis_prm_initialize(&mat->prm);
  monolis_com_initialize(&mat->com);
  monolis_mat_initialize(&mat->mat);
  monolis_mat_initialize(&mat->prec);
  monolis_com_input_comm_table(mat,
    mat->prm.com_top_dir_name,
    mat->prm.com_part_dir_name,
    mat->prm.com_file_name);
}

void monolis_initialize_entire(
  MONOLIS* mat)
{
  monolis_prm_initialize(&mat->prm);
  monolis_com_initialize(&mat->com);
  monolis_mat_initialize(&mat->mat);
  monolis_mat_initialize(&mat->prec);

  mat->com.comm = 0;
  mat->com.my_rank = 0;
  mat->com.comm_size = 1;
  mat->com.send_n_neib = 0;
  mat->com.recv_n_neib = 0;

  mat->com.send_neib_pe = monolis_alloc_I_1d(mat->com.send_neib_pe, 1);
  mat->com.send_index = monolis_alloc_I_1d(mat->com.send_index, 2);
  mat->com.send_item = monolis_alloc_I_1d(mat->com.send_item, 1);
  mat->com.recv_neib_pe = monolis_alloc_I_1d(mat->com.recv_neib_pe, 1);
  mat->com.recv_index = monolis_alloc_I_1d(mat->com.recv_index, 2);
  mat->com.recv_item = monolis_alloc_I_1d(mat->com.recv_item, 1);
}

void monolis_finalize(
  MONOLIS* mat)
{
  monolis_prm_initialize(&mat->prm);
  monolis_com_initialize(&mat->com);
  monolis_mat_initialize(&mat->mat);
  monolis_mat_initialize(&mat->prec);
}
