#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include "monolis_utils.h"
#include "monolis_def_struc_c.h"
#include "monolis_eigen_solver_c.h"

void monolis_eigen_standard_lanczos_R(
  MONOLIS* mat,
  int      n_get_eigen,
  double   ths,
  int      maxiter,
  double*  eigen_value,
  double** eigen_mode,
  bool*    is_Dirichlet_bc)
{
  int n = mat->mat.N;
  if(mat->com.comm_size > 1) n = mat->com.n_internal_vertex;
  int np = mat->mat.NP;
  int nz = mat->mat.CSR.index[np];
  int n_dof = mat->mat.NDOF;
  int recv_nitem = mat->com.recv_index[mat->com.recv_n_neib];
  int send_nitem = mat->com.send_index[mat->com.send_n_neib];
  int i, j;
  int* is_Dirichlet_bc_int;
  double* eigen_mode_tmp;

  eigen_mode_tmp = monolis_alloc_R_1d(eigen_mode_tmp, np*n_dof*n_get_eigen);
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
    mat->com.my_rank,
    mat->com.comm,
    mat->com.comm_size,
    mat->com.recv_n_neib,
    recv_nitem,
    mat->com.recv_neib_pe,
    mat->com.recv_index,
    mat->com.recv_item,
    mat->com.send_n_neib,
    send_nitem,
    mat->com.send_neib_pe,
    mat->com.send_index,
    mat->com.send_item,
    /* parameter */
    mat->prm.Iarray,
    mat->prm.Rarray,
    n_get_eigen,
    ths,
    maxiter,
    eigen_value,
    eigen_mode_tmp,
    is_Dirichlet_bc_int);

  for(i = 0; i < n_get_eigen; i++){
    for(j = 0; j < np*n_dof; j++){
      eigen_mode[i][j] = eigen_mode_tmp[np*n_dof*i + j];
    }
  }

  monolis_dealloc_R_1d(&eigen_mode_tmp);
  monolis_dealloc_I_1d(&is_Dirichlet_bc_int);
}

void monolis_eigen_inverted_standard_lanczos_R(
  MONOLIS* mat,
  int      n_get_eigen,
  double   ths,
  int      maxiter,
  double*  eigen_value,
  double** eigen_mode,
  bool*    is_Dirichlet_bc)
{
  int n = mat->mat.N;
  if(mat->com.comm_size > 1) n = mat->com.n_internal_vertex;
  int np = mat->mat.NP;
  int nz = mat->mat.CSR.index[np];
  int n_dof = mat->mat.NDOF;
  int recv_nitem = mat->com.recv_index[mat->com.recv_n_neib];
  int send_nitem = mat->com.send_index[mat->com.send_n_neib];
  int i, j;
  int* is_Dirichlet_bc_int;
  double* eigen_mode_tmp;

  eigen_mode_tmp = monolis_alloc_R_1d(eigen_mode_tmp, np*n_dof*n_get_eigen);
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
    mat->com.my_rank,
    mat->com.comm,
    mat->com.comm_size,
    mat->com.recv_n_neib,
    recv_nitem,
    mat->com.recv_neib_pe,
    mat->com.recv_index,
    mat->com.recv_item,
    mat->com.send_n_neib,
    send_nitem,
    mat->com.send_neib_pe,
    mat->com.send_index,
    mat->com.send_item,
    /* parameter */
    mat->prm.Iarray,
    mat->prm.Rarray,
    n_get_eigen,
    ths,
    maxiter,
    eigen_value,
    eigen_mode_tmp,
    is_Dirichlet_bc_int);

  for(i = 0; i < n_get_eigen; i++){
    for(j = 0; j < np*n_dof; j++){
      eigen_mode[i][j] = eigen_mode_tmp[np*n_dof*i + j];
    }
  }

  monolis_dealloc_R_1d(&eigen_mode_tmp);
  monolis_dealloc_I_1d(&is_Dirichlet_bc_int);
}
