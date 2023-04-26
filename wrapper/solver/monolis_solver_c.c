#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_def_struc_c.h"
#include "monolis_solver_c.h"

void monolis_solve_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  double*      b,
  double*      x)
{
  int n = mat->mat.N;
  if(com->comm_size > 1) n = com->n_internal_vertex;
  int np = mat->mat.NP;
  int nz = mat->mat.CSR.index[np];
  int n_dof = mat->mat.NDOF;
  int recv_nitem = com->recv_index[com->recv_n_neib];
  int send_nitem = com->send_index[com->send_n_neib];

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
    mat->prm.Rarray);
}

void monolis_solve_C(
  MONOLIS*         mat,
  MONOLIS_COM*     com,
  double _Complex* b,
  double _Complex* x)
{
  int n = mat->mat.N;
  if(com->comm_size > 1) n = com->n_internal_vertex;
  int np = mat->mat.NP;
  int nz = mat->mat.CSR.index[np];
  int n_dof = mat->mat.NDOF;
  int recv_nitem = com->recv_index[com->recv_n_neib];
  int send_nitem = com->send_index[com->send_n_neib];

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
    mat->prm.Rarray);
}
