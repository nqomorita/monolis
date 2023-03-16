#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_def_solver_c.h"
#include "monolis_def_struc_c.h"

void monolis_get_nonzero_pattern_by_nodal_graph_main(
  MONOLIS* mat,
  int      n_node,
  int      n_base,
  int      ndof,
  int      n_elem,
  int**    elem)
{
//  int i, j;
//  int nz, jS, jE;
//
//  mat->mat.N = nnode;
//  mat->mat.NP = nnode;
//  mat->mat.NDOF = ndof;
//
//  mat->mat.index = (int* )calloc(nnode+1, sizeof(int));
//  for(i = 1; i < nnode + 1; i++) {
//    mat->mat.index[i] = index[i] + i;
//  }
//
//  nz = mat->mat.index[nnode];
//  mat->mat.NZ = nz;
//  mat->mat.item = (int*)calloc(nz, sizeof(int));
//
//  for(i = 0; i < nnode; i++) {
//    jS = mat->mat.index[i];
//    jE = mat->mat.index[i+1];
//    mat->mat.item[jS] = i+1;
//    for(j = jS + 1; j < jE; j++){
//      mat->mat.item[j] = item[j-i-1] + 1;
//    }
//    monolis_qsort_int(
//      &(mat->mat.item[jS]),
//      1,
//      jE-jS);
//  }
//
//  mat->mat.indexR = (int*)calloc(nnode+1, sizeof(int));
//  mat->mat.itemR = (int*)calloc(nz, sizeof(int));
//  mat->mat.permR = (int*)calloc(nz, sizeof(int));
//  monolis_get_CRR_format(
//    nnode,
//    nnode,
//    nz,
//    mat->mat.index,
//    mat->mat.item,
//    mat->mat.indexR,
//    mat->mat.itemR,
//    mat->mat.permR);
}

void monolis_alloc_nonzero_pattern_mat_val_R()
{
//  mat->mat.A = (double*)calloc(ndof*ndof*nz, sizeof(double));
//  mat->mat.X = (double*)calloc(ndof*nnode, sizeof(double));
//  mat->mat.B = (double*)calloc(ndof*nnode, sizeof(double));
}

void monolis_alloc_nonzero_pattern_mat_val_C()
{

}
