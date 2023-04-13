/* monolis_spmat_handler.h */
#ifndef MONOLIS_SPMAT_HANDLER_H
#define MONOLIS_SPMAT_HANDLER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>

void monolis_set_scalar_to_sparse_matrix_R(
  MONOLIS* mat,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val);

void monolis_add_scalar_to_sparse_matrix_R(
  MONOLIS* mat,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val);

void monolis_get_scalar_from_sparse_matrix_R(
  MONOLIS* mat,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double*  val,
  bool*    is_find);

void monolis_add_matrix_to_sparse_matrix_R(
  MONOLIS* mat,
  int      n_base,
  int*     connectivity,
  double** val);

void monolis_add_matrix_to_sparse_matrix_offdiag_R(
  MONOLIS* mat,
  int      n_base1,
  int      n_base2,
  int*     connectivity1,
  int*     connectivity2,
  double** val);

void monolis_set_matrix_BCSR_R(
  MONOLIS* mat,
  int      n,
  int      np,
  int      n_dof,
  int      nz,
  double*  A,
  int*     index,
  int*     item);

void monolis_set_matrix_BCSR_mat_val_R(
  MONOLIS* mat,
  int      n_dof,
  int      nz,
  double*  A);

void monolis_set_Dirichlet_bc_R(
  MONOLIS* mat,
  double*  b,
  int      node_id,
  int      n_dof_bc,
  double   val);

void monolis_set_scalar_to_sparse_matrix_C(
  MONOLIS*       mat,
  int            i,
  int            j,
  int            submat_i,
  int            submat_j,
  double _Complex val);

void monolis_add_scalar_to_sparse_matrix_C(
  MONOLIS*       mat,
  int            i,
  int            j,
  int            submat_i,
  int            submat_j,
  double _Complex val);

void monolis_get_scalar_from_sparse_matrix_C(
  MONOLIS*        mat,
  int             i,
  int             j,
  int             submat_i,
  int             submat_j,
  double _Complex* val,
  bool*           is_find);

void monolis_add_matrix_to_sparse_matrix_C(
  MONOLIS*         mat,
  int              n_base,
  int*             connectivity,
  double _Complex** val);

void monolis_add_matrix_to_sparse_matrix_offdiag_C(
  MONOLIS*         mat,
  int              n_base1,
  int              n_base2,
  int*             connectivity1,
  int*             connectivity2,
  double _Complex** val);

void monolis_set_matrix_BCSR_C(
  MONOLIS*        mat,
  int             n,
  int             np,
  int             n_dof,
  int             nz,
  double _Complex* A,
  int*            index,
  int*            item);

void monolis_set_matrix_BCSR_mat_val_C(
  MONOLIS*        mat,
  int             n_dof,
  int             nz,
  double _Complex* A);

void monolis_set_Dirichlet_bc_C(
  MONOLIS*        mat,
  double _Complex* b,
  int             node_id,
  int             n_dof_bc,
  double _Complex  val);

#ifdef __cplusplus
}
#endif

#endif
