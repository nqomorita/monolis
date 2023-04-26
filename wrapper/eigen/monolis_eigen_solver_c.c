#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis_utils.h"
#include "monolis_def_struc_c.h"
#include "monolis_eigen_solver_c.h"

void monolis_eigen_standard_lanczos_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  int*         n_get_eigen,
  double       ths,
  int          maxiter,
  double*      eigen_value,
  double**     eigen_mode,
  bool*        is_Dirichlet_bc)
{
  int n = mat->mat.N;
  if(com->comm_size > 1) n = com->n_internal_vertex;
  int np = mat->mat.NP;
  int nz = mat->mat.CSR.index[np];
  int n_dof = mat->mat.NDOF;
  int recv_nitem = com->recv_index[com->recv_n_neib];
  int send_nitem = com->send_index[com->send_n_neib];
  int i, j;
  int* is_Dirichlet_bc_int;
  double* eigen_mode_tmp;

  eigen_mode_tmp = monolis_alloc_R_1d(eigen_mode_tmp, np*n_dof*(*n_get_eigen));
  is_Dirichlet_bc_int = monolis_alloc_I_1d(is_Dirichlet_bc_int, np*n_dof);

  for(i = 0; i < np*n_dof; i++){
    if (is_Dirichlet_bc[i]) {
      is_Dirichlet_bc_int[i] = 1;
    } else {
      is_Dirichlet_bc_int[i] = 0;
    }
  }

  monolis_eigen_standard_lanczos_R_c_main(
    /* mat */
    n,
    np,
    nz,
    n_dof,
    mat->mat.R.A,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    /* comm */
    com->my_rank,
    com->comm,
    com->comm_size,
    com->recv_n_neib,
    recv_nitem,
    com->recv_neib_pe,
    com->recv_index,
    com->recv_item,
    com->send_n_neib,
    send_nitem,
    com->send_neib_pe,
    com->send_index,
    com->send_item,
    /* parameter */
    mat->prm.Iarray,
    mat->prm.Rarray,
    n_get_eigen,
    ths,
    maxiter,
    eigen_value,
    eigen_mode_tmp,
    is_Dirichlet_bc_int);

  for(i = 0; i < *n_get_eigen; i++){
    for(j = 0; j < np*n_dof; j++){
      eigen_mode[i][j] = eigen_mode_tmp[np*n_dof*i + j];
    }
  }

  monolis_dealloc_R_1d(&eigen_mode_tmp);
  monolis_dealloc_I_1d(&is_Dirichlet_bc_int);
}

void monolis_eigen_inverted_standard_lanczos_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  int*         n_get_eigen,
  double       ths,
  int          maxiter,
  double*      eigen_value,
  double**     eigen_mode,
  bool*        is_Dirichlet_bc)
{
  int n = mat->mat.N;
  if(com->comm_size > 1) n = com->n_internal_vertex;
  int np = mat->mat.NP;
  int nz = mat->mat.CSR.index[np];
  int n_dof = mat->mat.NDOF;
  int recv_nitem = com->recv_index[com->recv_n_neib];
  int send_nitem = com->send_index[com->send_n_neib];
  int i, j;
  int* is_Dirichlet_bc_int;
  double* eigen_mode_tmp;

  eigen_mode_tmp = monolis_alloc_R_1d(eigen_mode_tmp, np*n_dof*(*n_get_eigen));
  is_Dirichlet_bc_int = monolis_alloc_I_1d(is_Dirichlet_bc_int, np*n_dof);

  for(i = 0; i < np*n_dof; i++){
    if (is_Dirichlet_bc[i]) {
      is_Dirichlet_bc_int[i] = 1;
    } else {
      is_Dirichlet_bc_int[i] = 0;
    }
  }

  monolis_eigen_inverted_standard_lanczos_R_c_main(
    /* mat */
    n,
    np,
    nz,
    n_dof,
    mat->mat.R.A,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    /* comm */
    com->my_rank,
    com->comm,
    com->comm_size,
    com->recv_n_neib,
    recv_nitem,
    com->recv_neib_pe,
    com->recv_index,
    com->recv_item,
    com->send_n_neib,
    send_nitem,
    com->send_neib_pe,
    com->send_index,
    com->send_item,
    /* parameter */
    mat->prm.Iarray,
    mat->prm.Rarray,
    n_get_eigen,
    ths,
    maxiter,
    eigen_value,
    eigen_mode_tmp,
    is_Dirichlet_bc_int);

  for(i = 0; i < *n_get_eigen; i++){
    for(j = 0; j < np*n_dof; j++){
      eigen_mode[i][j] = eigen_mode_tmp[np*n_dof*i + j];
    }
  }

  monolis_dealloc_R_1d(&eigen_mode_tmp);
  monolis_dealloc_I_1d(&is_Dirichlet_bc_int);
}
