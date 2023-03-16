/* monolis_spmat_handler.h */
#ifndef MONOLIS_SPMAT_HANDLER_H
#define MONOLIS_SPMAT_HANDLER_H

#ifdef __cplusplus
extern "C" {
#endif

void monolis_set_scalar_to_sparse_matrix_R_c_main(
  int      n_node,
  int      nz,
  int      n_dof,
  int*     index,
  int*     item,
  double*  A,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val);

void monolis_add_scalar_to_sparse_matrix_R_c_main(
  int      n_node,
  int      nz,
  int      n_dof,
  int*     index,
  int*     item,
  double*  A,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val);

void monolis_get_scalar_from_sparse_matrix_R_c_main(
  int      n_node,
  int      nz,
  int      n_dof,
  int*     index,
  int*     item,
  double*  A,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double*  val,
  int*     is_find);

void monolis_add_matrix_to_sparse_matrix_main_R_c_main(
  int     n_node,
  int     nz,
  int     n_dof,
  int*    index,
  int*    item,
  double* A,
  int     n_base1,
  int     n_base2,
  int*    connectivity1,
  int*    connectivity2,
  double* val);

void monolis_set_Dirichlet_bc_R_c_main(
  int      n_node,
  int      nz,
  int      n_dof,
  int*     index,
  int*     item,
  int*     indexR,
  int*     itemR,
  int*     permR,
  double*  A,
  double*  b,
  int      node_id,
  int      n_dof_bc,
  double   val);

void monolis_set_scalar_to_sparse_matrix_C_c_main(
  int             n_node,
  int             nz,
  int             n_dof,
  int*            index,
  int*            item,
  double complex* A,
  int             i,
  int             j,
  int             submat_i,
  int             submat_j,
  double complex  val);

void monolis_add_scalar_to_sparse_matrix_C_c_main(
  int             n_node,
  int             nz,
  int             n_dof,
  int*            index,
  int*            item,
  double complex* A,
  int             i,
  int             j,
  int             submat_i,
  int             submat_j,
  double complex  val);

void monolis_get_scalar_from_sparse_matrix_C_c_main(
  int             n_node,
  int             nz,
  int             n_dof,
  int*            index,
  int*            item,
  double complex* A,
  int             i,
  int             j,
  int             submat_i,
  int             submat_j,
  double complex* val,
  int*            is_find);

void monolis_add_matrix_to_sparse_matrix_main_C_c_main(
  int              n_node,
  int              nz,
  int              n_dof,
  int*             index,
  int*             item,
  double complex*  A,
  int              n_base1,
  int              n_base2,
  int*             connectivity1,
  int*             connectivity2,
  double complex** val);

void monolis_set_Dirichlet_bc_C_c_main(
  int             n_node,
  int             nz,
  int             n_dof,
  int*            index,
  int*            item,
  int*            indexR,
  int*            itemR,
  int*            permR,
  double complex* A,
  double complex* b,
  int             node_id,
  int             n_dof_bc,
  double          val);

#ifdef __cplusplus
}
#endif

#endif
