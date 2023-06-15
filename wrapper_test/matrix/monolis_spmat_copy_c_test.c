#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_spmat_copy_c_test_R()
{
  MONOLIS mat_in;
  MONOLIS mat_out;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int i, nz;
  int** elem;

  n_node = 4;
  n_base = 2;
  n_dof  = 2;
  n_elem = 3;

  elem = monolis_alloc_I_2d(elem, n_node, n_base);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;

  monolis_get_nonzero_pattern_by_simple_mesh_R(&mat_in, n_node, n_base, n_dof, n_elem, elem);

  nz = mat_in.mat.CSR.index[n_node];
  for (int i = 0; i < nz; ++i) {
    mat_in.mat.R.A[i] = 1.0;
  }

  for (int i = 0; i < n_node*n_dof; ++i) {
    mat_in.mat.R.B[i] = 2.0;
  }

  for (int i = 0; i < n_node*n_dof; ++i) {
    mat_in.mat.R.X[i] = 3.0;
  }

  monolis_copy_mat_R(&mat_in, &mat_out);

  //monolis_test_check_eq_R1("monolis_copy_mat_R", mat_out.mat.R.A[12], 3.0);



  //monolis_clear_mat_value_R(&mat_out);

  for (int i = 0; i < nz; ++i) {
    mat_in.mat.R.A[i] = 0.0;
  }

  for (int i = 0; i < n_node*n_dof; ++i) {
    mat_in.mat.R.B[i] = 0.0;
  }

  for (int i = 0; i < n_node*n_dof; ++i) {
    mat_in.mat.R.X[i] = 0.0;
  }
}

void monolis_spmat_copy_c_test_C()
{

}

void monolis_spmat_copy_c_test()
{
  monolis_std_log_string("monolis_spmat_copy_c_test");
  monolis_spmat_copy_c_test_R();
  monolis_spmat_copy_c_test_C();
}
