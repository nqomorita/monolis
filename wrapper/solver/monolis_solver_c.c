#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "monolis_utils.h"
#include "monolis_def_struc_c.h"
#include "monolis_solver_c.h"

void monolis_solve_R(
  MONOLIS* mat,
  double*  b,
  double*  x)
{
  int n = mat->mat.N;
  if(mat->com.comm_size > 1) n = mat->com.n_internal_vertex;
  int np = mat->mat.NP;
  int nz = mat->mat.CSR.index[np];
  int n_dof = mat->mat.NDOF;
  int recv_nitem = mat->com.recv_index[mat->com.recv_n_neib];
  int send_nitem = mat->com.send_index[mat->com.send_n_neib];

  monolis_solve_R_c_main(
    /* mat */
    n,
    np,
    nz,
    n_dof,
    mat->mat.R.A,
    x,
    b,
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
    mat->prm.Rarray);
}

void monolis_solve_C(
  MONOLIS* mat,
  double complex*  b,
  double complex*  x)
{
  int n = mat->mat.N;
  if(mat->com.comm_size > 1) n = mat->com.n_internal_vertex;
  int np = mat->mat.NP;
  int nz = mat->mat.CSR.index[np];
  int n_dof = mat->mat.NDOF;
  int recv_nitem = mat->com.recv_index[mat->com.recv_n_neib];
  int send_nitem = mat->com.send_index[mat->com.send_n_neib];

  monolis_solve_C_c_main(
    /* mat */
    n,
    np,
    nz,
    n_dof,
    mat->mat.C.A,
    x,
    b,
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
    mat->prm.Rarray);
}
