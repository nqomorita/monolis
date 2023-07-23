#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "gedatsu.h"
#include "monolis_spmat_nzpattern_c.h"
#include "monolis_spmat_nzpattern_util_c.h"
#include "monolis_def_struc_c.h"

void monolis_get_nonzero_pattern_by_simple_mesh_R(
  MONOLIS* mat,
  int      n_node,
  int      n_base,
  int      n_dof,
  int      n_elem,
  int**    elem)
{
  int* conn_index;
  int* conn_item;
  int* index;
  int* item;

  conn_index = monolis_alloc_I_1d(conn_index, n_elem + 1);
  conn_item = monolis_alloc_I_1d(conn_item,  n_elem*n_base);

  gedatsu_convert_simple_mesh_to_connectivity_graph(
    n_elem,
    n_base,
    elem,
    conn_index,
    conn_item);

  gedatsu_convert_connectivity_graph_to_nodal_graph(
    n_node,
    n_elem,
    conn_index,
    conn_item,
    &index,
    &item);

  monolis_get_nonzero_pattern_by_nodal_graph_main(
    &mat->mat,
    n_node,
    n_dof,
    index,
    item);

  monolis_alloc_nonzero_pattern_mat_val_R(
    &mat->mat);

  monolis_dealloc_I_1d(&conn_index);
  monolis_dealloc_I_1d(&conn_item);
  monolis_dealloc_I_1d(&index);
  monolis_dealloc_I_1d(&item);
}

void monolis_get_nonzero_pattern_by_connectivity_R(
  MONOLIS* mat,
  int      n_node,
  int      n_dof,
  int      n_elem,
  int*     conn_index,
  int*     conn_item)
{
  int* index;
  int* item;

  gedatsu_convert_connectivity_graph_to_nodal_graph(
    n_node,
    n_elem,
    conn_index,
    conn_item,
    &index,
    &item);

  monolis_get_nonzero_pattern_by_nodal_graph_main(
    &mat->mat,
    n_node,
    n_dof,
    index,
    item);

  monolis_alloc_nonzero_pattern_mat_val_R(
    &mat->mat);

  monolis_dealloc_I_1d(&conn_index);
  monolis_dealloc_I_1d(&conn_item);
  monolis_dealloc_I_1d(&index);
  monolis_dealloc_I_1d(&item);
}

void monolis_get_nonzero_pattern_by_nodal_graph_R(
  MONOLIS* mat,
  int      n_node,
  int      n_dof,
  int*     index,
  int*     item)
{
  monolis_get_nonzero_pattern_by_nodal_graph_main(
    &mat->mat,
    n_node,
    n_dof,
    index,
    item);

  monolis_alloc_nonzero_pattern_mat_val_R(
    &mat->mat);
}

void monolis_get_nonzero_pattern_by_simple_mesh_C(
  MONOLIS* mat,
  int      n_node,
  int      n_base,
  int      n_dof,
  int      n_elem,
  int**    elem)
{
  int* conn_index;
  int* conn_item;
  int* index;
  int* item;

  conn_index = monolis_alloc_I_1d(conn_index, n_elem + 1);
  conn_item = monolis_alloc_I_1d(conn_item,  n_elem*n_base);

  gedatsu_convert_simple_mesh_to_connectivity_graph(
    n_elem,
    n_base,
    elem,
    conn_index,
    conn_item);

  gedatsu_convert_connectivity_graph_to_nodal_graph(
    n_node,
    n_elem,
    conn_index,
    conn_item,
    &index,
    &item);

  monolis_get_nonzero_pattern_by_nodal_graph_main(
    &mat->mat,
    n_node,
    n_dof,
    index,
    item);

  monolis_alloc_nonzero_pattern_mat_val_C(
    &mat->mat);

  monolis_dealloc_I_1d(&conn_index);
  monolis_dealloc_I_1d(&conn_item);
  monolis_dealloc_I_1d(&index);
  monolis_dealloc_I_1d(&item);
}

void monolis_get_nonzero_pattern_by_connectivity_C(
  MONOLIS* mat,
  int      n_node,
  int      n_dof,
  int      n_elem,
  int*     conn_index,
  int*     conn_item)
{
  int* index;
  int* item;

  gedatsu_convert_connectivity_graph_to_nodal_graph(
    n_node,
    n_elem,
    conn_index,
    conn_item,
    &index,
    &item);

  monolis_get_nonzero_pattern_by_nodal_graph_main(
    &mat->mat,
    n_node,
    n_dof,
    index,
    item);

  monolis_alloc_nonzero_pattern_mat_val_C(
    &mat->mat);

  monolis_dealloc_I_1d(&conn_index);
  monolis_dealloc_I_1d(&conn_item);
  monolis_dealloc_I_1d(&index);
  monolis_dealloc_I_1d(&item);
}

void monolis_get_nonzero_pattern_by_nodal_graph_C(
  MONOLIS* mat,
  int      n_node,
  int      n_dof,
  int*     index,
  int*     item)
{
  monolis_get_nonzero_pattern_by_nodal_graph_main(
    &mat->mat,
    n_node,
    n_dof,
    index,
    item);

  monolis_alloc_nonzero_pattern_mat_val_C(
    &mat->mat);
}
