#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_matvec_c.h"
#include "monolis_def_struc_c.h"

void monolis_matvec_product_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  double*      x,
  double*      y)
{
  int n = mat->mat.N;
  if(com->comm_size > 1 ) n = com->n_internal_vertex;
  int np = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[mat->mat.NP];
  int recv_nitem = com->recv_index[com->recv_n_neib];
  int send_nitem = com->send_index[com->send_n_neib];

  monolis_matvec_product_R_c_main(
    /* mat */
    n,
    np,
    nz,
    ndof,
    mat->mat.R.A,
    x,
    y,
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
    com->send_item);
}

void monolis_matvec_product_C(
  MONOLIS*        mat,
  MONOLIS_COM*    com,
  double _Complex* x,
  double _Complex* y)
{
  int n = mat->mat.N;
  if(com->comm_size > 1 ) n = com->n_internal_vertex;
  int np = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[mat->mat.NP];
  int recv_nitem = com->recv_index[com->recv_n_neib];
  int send_nitem = com->send_index[com->send_n_neib];

  monolis_matvec_product_C_c_main(
    /* mat */
    n,
    np,
    nz,
    ndof,
    mat->mat.C.A,
    x,
    y,
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
    com->send_item);
}
