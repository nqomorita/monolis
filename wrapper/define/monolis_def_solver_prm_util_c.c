#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_def_solver_prm_util_c.h"
#include "monolis_def_solver_prm_c.h"


void monolis_set_input_top_directory_name(
  MONOLIS*    mat,
  const char* param)
{
  strcpy(mat->prm.com_top_dir_name, param);
}

void monolis_set_input_part_directory_name(
  MONOLIS*    mat,
  const char* param)
{
  strcpy(mat->prm.com_part_dir_name, param);
}

void monolis_set_input_file_name(
  MONOLIS*    mat,
  const char* param)
{
  strcpy(mat->prm.com_file_name, param);
}
















void monolis_set_communicator(
  MONOLIS* mat,
  int      param)
{
  monolis_com_set_communicator(&mat->com, param);
}

void monolis_set_my_rank(
  MONOLIS* mat,
  int      param)
{
  monolis_com_set_my_rank(&mat->com, param);
}

void monolis_set_comm_size(
  MONOLIS* mat,
  int      param)
{
  monolis_com_set_comm_size(&mat->com, param);
}

void monolis_set_n_internal_vertex(
  MONOLIS* mat,
  int      param)
{
  monolis_com_set_n_internal_vertex(&mat->com, param);
}

void monolis_get_n_internal_vertex(
  MONOLIS* mat,
  int*     param)
{
  monolis_com_get_n_internal_vertex(&mat->com, param);
}

void monolis_get_internal_simple_mesh_list(
  MONOLIS* mat,
  int      n_elem,
  int      n_base,
  int**    elem,
  int*     list)
{
  int i, j;
  int my_rank;
  int n_internal_vertex;
  int* domain_id;

  monolis_get_n_internal_vertex(mat, &n_internal_vertex);

  //my_rank = monolis_get_local_my_rank(&mat);

  monolis_alloc_I_1d(domain_id, mat->mat.NP);

  for (i = 0; i < mat->mat.NP; ++i) {
    domain_id[i] = my_rank;
  }

  //monolis_update_I(mat->com, 1, domain_id);

  for (i = 0; i < mat->mat.NP; ++i) {
    list[i] = 0;
  }

  for (i = 0; i < n_elem; ++i) {
    for (j = 0; i < n_base; ++i) {

    }
  }


/*
    do i = 1, n_elem
      do j = 1, n_base
        id(j) = domain_id(elem(j,i))
      enddo
      in = minval(id)
      if(in == my_rank) list(i) = 1
    enddo
*/
}

void monolis_get_internal_connectivity_list(
  MONOLIS* mat,
  int      n_elem,
  int*     index,
  int*     item,
  int*     list)
{
/*
    call monolis_get_n_internal_vertex(monolis, n_internal_vertex)

    my_rank = monolis_get_local_my_rank(monolis)

    call monolis_alloc_I_1d(domain_id, monolis%MAT%NP)

    do i = 1, n_internal_vertex
      domain_id(i) = my_rank
    enddo

    call monolis_mpi_update_I(monolis%COM, 1, domain_id)

    list = 0
    do i = 1, n_elem
      jS = index(i) + 1
      jE = index(i + 1)
      call monolis_alloc_I_1d(id, jE - jS + 1)
      do j = jS, jE
        id(j - jS + 1) = domain_id(item(j))
      enddo
      call monolis_dealloc_I_1d(id)
      in = minval(id)
      if(in == my_rank) list(i) = 1
    enddo
*/
}

void monolis_set_method(
  MONOLIS* mat,
  int      param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_METHOD] = param;
}

void monolis_set_precond(
  MONOLIS* mat,
  int      param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_PRECOND] = param;
}

void monolis_set_maxiter(
  MONOLIS* mat,
  int      param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_MAX_ITER] = param;
}

void monolis_get_converge_iter(
  MONOLIS* mat,
  int*     param)
{
  *param = mat->prm.Iarray[MONOLIS_PRM_I_CUR_ITER];
}

void monolis_get_error_tag(
  MONOLIS* mat,
  int*     param)
{
  *param = mat->prm.Iarray[MONOLIS_PRM_I_IERR];
}

void monolis_set_init_x(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_INIT_X] = (int)param;
}

void monolis_set_sym_matrix(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_SYM_MATRIX] = param;
}

void monolis_set_debug(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_DEBUG] = param;
}

void monolis_set_performance_measurement(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_MEASUREMENT] = param;
}

void monolis_set_check_diag(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_CHECK_DIAG] = param;
}

void monolis_set_prec_stored(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_PREC_STORED] = param;
}

void monolis_show_iterlog(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_SHOW_ITERLOG] = param;
}

void monolis_show_timelog(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_SHOW_TIME] = param;
}

void monolis_show_summary(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_SHOW_SUMMARY] = param;
}

void monolis_show_timelog_statistics(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_SHOW_TIME_STATISTICS] = param;
}

void monolis_set_tolerance(
  MONOLIS* mat,
  double   val)
{
  mat->prm.Rarray[MONOLIS_PRM_R_TOL] = val;
}

void monolis_get_converge_residual(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_PRM_R_CUR_RESID];
}

void monolis_get_time_solver(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_SOL];
}

void monolis_get_time_preparing(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_PREP];
}

void monolis_get_time_spmv(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_SPMV];
}

void monolis_get_time_inner_product(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_DOTP];
}

void monolis_get_time_precondition(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_PREC];
}

void monolis_get_time_comm_inner_product(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_COMM_DOTP];
}

void monolis_get_time_comm_spmv(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_COMM_SPMV];
}
