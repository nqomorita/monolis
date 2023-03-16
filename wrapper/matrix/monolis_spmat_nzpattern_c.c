#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "gedatsu.h"
#include "monolis_spmat_nzpattern_c.h"
#include "monolis_def_struc_c.h"

void monolis_get_nonzero_pattern_by_simple_mesh_R(
  MONOLIS* mat,
  int      n_node,
  int      n_base,
  int      ndof,
  int      n_elem,
  int**    elem)
{
  int* conn_index;
  int* con;
  int* index;
  int* item;

  //conn_index = (idx_t*)calloc(n_elem+1, sizeof(idx_t));
  //con = (idx_t*)calloc(n_elem*n_base, sizeof(idx_t));

  gedatsu_convert_simple_mesh_to_connectivity_graph(
    n_elem,
    n_base,
    elem,
    conn_index,
    con);

  gedatsu_convert_connectivity_graph_to_nodal_graph(
    n_node,
    n_elem,
    conn_index,
    con,
    &index,
    &item);

//  monolis_get_nonzero_pattern_by_nodal_graph_main(
//    mat,
//    n_node,
//    ndof,
//    index,
//    item);

//  monolis_alloc_nonzero_pattern_mat_val_R();

  free(conn_index);
  free(con);
  free(index);
  free(item);
}

void monolis_get_nonzero_pattern_by_connectivity_R(
  MONOLIS* mat,
  int      n_node,
  int      n_base,
  int      ndof,
  int      n_elem,
  int**    elem)
{
  int* conn_index;
  int* con;
  int* index;
  int* item;

  gedatsu_convert_connectivity_graph_to_nodal_graph(
    n_node,
    n_elem,
    conn_index,
    con,
    &index,
    &item);

//  monolis_get_nonzero_pattern_by_nodal_graph_main(
//    mat,
//    n_node,
//    ndof,
//    index,
//    item);

//  monolis_alloc_nonzero_pattern_mat_val_R();
}

void monolis_get_nonzero_pattern_by_nodal_graph_R(
  MONOLIS* mat,
  int      n_node,
  int      n_base,
  int      ndof,
  int      n_elem,
  int**    elem)
{
//  monolis_get_nonzero_pattern_by_nodal_graph_main(
//    mat,
//    n_node,
//    ndof,
//    index,
//    item);

//  monolis_alloc_nonzero_pattern_mat_val_R();
}
