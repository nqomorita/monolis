#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "monolis_utils.h"
#include "monolis_matvec_c.h"
#include "monolis_def_struc_c.h"

void monolis_matvec_product_R(
  MONOLIS* mat,
  double*  x,
  double*  y)
{
  int n = mat->mat.N;
  if(mat->com.comm_size > 1 ) n = mat->com.n_internal_vertex;
  int np = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[mat->mat.NP + 1];
  int recv_nitem = mat->com.recv_index[mat->com.recv_n_neib];
  int send_nitem = mat->com.send_index[mat->com.send_n_neib];

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
    mat->com.send_item);
}

void monolis_matvec_product_C(
  MONOLIS*        mat,
  double complex* x,
  double complex* y)
{
  int n = mat->mat.N;
  if(mat->com.comm_size > 1 ) n = mat->com.n_internal_vertex;
  int np = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[mat->mat.NP + 1];
  int recv_nitem = mat->com.recv_index[mat->com.recv_n_neib];
  int send_nitem = mat->com.send_index[mat->com.send_n_neib];

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
    mat->com.send_item);
}
