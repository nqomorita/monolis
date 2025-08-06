!> 疎行列操作関数群（メイン関数）
module mod_monolis_spmat_handler_util_wrap
  use mod_monolis_utils
  use mod_monolis_spmat_handler
  use mod_monolis_spmat_handler_util
  use iso_c_binding

  implicit none

contains

  subroutine monolis_set_scalar_to_sparse_matrix_R_c(N, NZ, NZA, n_dof_index, n_dof_index2, &
    index, item, A, i, j, sub_i, sub_j, val) &
    bind(c, name = "monolis_set_scalar_to_sparse_matrix_R_c_main")
    implicit none
    integer(kint_c), intent(in), value :: N, NZ, NZA, i, j, sub_i, sub_j
    integer(kint_c), intent(in) :: index(N + 1)
    integer(kint_c), intent(in) :: item(NZ)
    integer(kint_c), intent(in) :: n_dof_index(N + 1)
    integer(kint_c), intent(in) :: n_dof_index2(NZ)
    real(c_double) :: A(NZA)
    real(c_double), value :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_set_scalar_to_sparse_matrix_main_R(index, item, A, n_dof_index, n_dof_index2, &
      i_t, j_t, sub_i_t, sub_j_t, val)
  end subroutine monolis_set_scalar_to_sparse_matrix_R_c

  subroutine monolis_add_scalar_to_sparse_matrix_R_c(N, NZ, NZA, n_dof_index, n_dof_index2, &
    index, item, A, i, j, sub_i, sub_j, val) &
    bind(c, name = "monolis_add_scalar_to_sparse_matrix_R_c_main")
    implicit none
    integer(kint_c), intent(in), value :: N, NZ, NZA, i, j, sub_i, sub_j
    integer(kint_c), intent(in) :: index(N + 1)
    integer(kint_c), intent(in) :: item(NZ)
    integer(kint_c), intent(in) :: n_dof_index(N + 1)
    integer(kint_c), intent(in) :: n_dof_index2(NZ)
    real(c_double) :: A(NZA)
    real(c_double), value :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_add_scalar_to_sparse_matrix_main_R(index, item, A, n_dof_index, n_dof_index2, &
      i_t, j_t, sub_i_t, sub_j_t, val)
  end subroutine monolis_add_scalar_to_sparse_matrix_R_c

  subroutine monolis_get_scalar_from_sparse_matrix_R_c(N, NZ, NZA, n_dof_index, n_dof_index2, &
    index, item, A, i, j, sub_i, sub_j, val, is_find) &
    bind(c, name = "monolis_get_scalar_from_sparse_matrix_R_c_main")
    implicit none
    integer(kint_c), intent(in), value :: N, NZ, NZA, i, j, sub_i, sub_j
    integer(kint_c), intent(in) :: index(N + 1)
    integer(kint_c), intent(in) :: item(NZ)
    integer(kint_c), intent(in) :: n_dof_index(N + 1)
    integer(kint_c), intent(in) :: n_dof_index2(NZ)
    real(c_double) :: A(NZA)
    integer(kint_c) :: is_find
    real(c_double) :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t
    logical :: is_find_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_get_scalar_from_sparse_matrix_main_R(index, item, A, n_dof_index, n_dof_index2, &
      i_t, j_t, sub_i_t, sub_j_t, val, is_find_t)
    is_find = 0
    if(is_find_t) is_find = 1
  end subroutine monolis_get_scalar_from_sparse_matrix_R_c

  subroutine monolis_add_matrix_to_sparse_matrix_main_R_c(N, NZ, NZA, NZM1, NZM2, n_dof_index, n_dof_index2, &
    index, item, A, n_base1, n_base2, conn1, conn2, mat) &
    bind(c, name = "monolis_add_matrix_to_sparse_matrix_main_R_c_main")
    implicit none
    integer(kint_c), intent(in), value :: N, NZ, NZA, NZM1, NZM2, n_base1, n_base2
    integer(kint_c), intent(in) :: index(N + 1)
    integer(kint_c), intent(in) :: item(NZ)
    integer(kint_c), intent(in) :: conn1(n_base1)
    integer(kint_c), intent(in) :: conn2(n_base2)
    integer(kint_c), intent(in) :: n_dof_index(N + 1)
    integer(kint_c), intent(in) :: n_dof_index2(NZ)
    real(c_double) :: A(NZA)
    real(c_double) :: mat(NZM1*NZM2)
    integer(kint) :: conn1t(n_base1), conn2t(n_base2), i, j
    real(kdouble) :: mat_t(NZM1,NZM2)

    conn1t = conn1 + 1
    conn2t = conn2 + 1

    do i = 1, NZM1
      do j = 1, NZM2
        mat_t(i,j) = mat(NZM2*(i - 1) + j)
      enddo
    enddo

    call monolis_add_matrix_to_sparse_matrix_main_R(index, item, A, n_base1, n_base2, &
      n_dof_index, n_dof_index2, conn1t, conn2t, mat_t)
  end subroutine monolis_add_matrix_to_sparse_matrix_main_R_c

  subroutine monolis_set_Dirichlet_bc_R_c(N, NZ, NZA, NZB, n_dof_index, n_dof_index2, &
    index, item, indexR, itemR, permR, A, B, node_id, ndof_bc, val) &
    bind(c, name = "monolis_set_Dirichlet_bc_R_c_main")
    implicit none
    integer(kint_c), intent(in), value :: N, NZ, NZA, NZB, node_id, ndof_bc
    integer(kint_c), intent(in) :: index(N + 1)
    integer(kint_c), intent(in) :: item(NZ)
    integer(kint_c), intent(in) :: indexR(N + 1)
    integer(kint_c), intent(in) :: itemR(NZ)
    integer(kint_c), intent(in) :: permR(NZ)
    real(c_double), intent(in), value :: val
    integer(kint_c), intent(in) :: n_dof_index(N + 1)
    integer(kint_c), intent(in) :: n_dof_index2(NZ)
    real(c_double) :: A(NZA)
    real(c_double) :: B(NZB)
    integer(kint) :: nid_t, ndof_bc_t

    nid_t = node_id + 1
    ndof_bc_t = ndof_bc + 1
    call monolis_set_Dirichlet_bc_main_R(index, item, A, B, indexR, itemR, permR, &
      & n_dof_index, n_dof_index2, nid_t, ndof_bc_t, val)
  end subroutine monolis_set_Dirichlet_bc_R_c

  subroutine monolis_set_scalar_to_sparse_matrix_C_c(N, NZ, NZA, n_dof_index, n_dof_index2, &
    index, item, A, i, j, sub_i, sub_j, val) &
    bind(c, name = "monolis_set_scalar_to_sparse_matrix_C_c_main")
    implicit none
    integer(kint_c), intent(in), value :: N, NZ, NZA, i, j, sub_i, sub_j
    integer(kint_c), intent(in) :: index(N + 1)
    integer(kint_c), intent(in) :: item(NZ)
    integer(kint_c), intent(in) :: n_dof_index(N + 1)
    integer(kint_c), intent(in) :: n_dof_index2(NZ)
    complex(c_double) :: A(NZA)
    complex(c_double), value :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_set_scalar_to_sparse_matrix_main_C(index, item, A, n_dof_index, n_dof_index2, &
      i_t, j_t, sub_i_t, sub_j_t, val)
  end subroutine monolis_set_scalar_to_sparse_matrix_C_c

  subroutine monolis_add_scalar_to_sparse_matrix_C_c(N, NZ, NZA, n_dof_index, n_dof_index2, &
    index, item, A, i, j, sub_i, sub_j, val) &
    bind(c, name = "monolis_add_scalar_to_sparse_matrix_C_c_main")
    implicit none
    integer(kint_c), intent(in), value :: N, NZ, NZA, i, j, sub_i, sub_j
    integer(kint_c), intent(in) :: index(N + 1)
    integer(kint_c), intent(in) :: item(NZ)
    integer(kint_c), intent(in) :: n_dof_index(N + 1)
    integer(kint_c), intent(in) :: n_dof_index2(NZ)
    complex(c_double) :: A(NZA)
    complex(c_double), value :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_add_scalar_to_sparse_matrix_main_C(index, item, A, n_dof_index, n_dof_index2, &
      i_t, j_t, sub_i_t, sub_j_t, val)
  end subroutine monolis_add_scalar_to_sparse_matrix_C_c

  subroutine monolis_get_scalar_from_sparse_matrix_C_c(N, NZ, NZA, n_dof_index, n_dof_index2, &
    index, item, A, i, j, sub_i, sub_j, val, is_find) &
    bind(c, name = "monolis_get_scalar_from_sparse_matrix_C_c_main")
    implicit none
    integer(kint_c), intent(in), value :: N, NZ, NZA, i, j, sub_i, sub_j
    integer(kint_c), intent(in) :: index(N + 1)
    integer(kint_c), intent(in) :: item(NZ)
    integer(kint_c) :: is_find
    integer(kint_c), intent(in) :: n_dof_index(N + 1)
    integer(kint_c), intent(in) :: n_dof_index2(NZ)
    complex(c_double) :: A(NZA)
    complex(c_double) :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t
    logical :: is_find_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_get_scalar_from_sparse_matrix_main_C(index, item, A, n_dof_index, n_dof_index2, &
      i_t, j_t, sub_i_t, sub_j_t, val, is_find_t)
    is_find = 0
    if(is_find_t) is_find = 1
  end subroutine monolis_get_scalar_from_sparse_matrix_C_c

  subroutine monolis_add_matrix_to_sparse_matrix_main_C_c(N, NZ, NZA, NZM1, NZM2, n_dof_index, n_dof_index2, &
    index, item, A, n_base1, n_base2, conn1, conn2, mat) &
    bind(c, name = "monolis_add_matrix_to_sparse_matrix_main_C_c_main")
    implicit none
    integer(kint_c), intent(in), value :: N, NZ, NZA, NZM1, NZM2, n_base1, n_base2
    integer(kint_c), intent(in) :: index(N + 1)
    integer(kint_c), intent(in) :: item(NZ)
    integer(kint_c), intent(in) :: conn1(n_base1)
    integer(kint_c), intent(in) :: conn2(n_base2)
    integer(kint_c), intent(in) :: n_dof_index(N + 1)
    integer(kint_c), intent(in) :: n_dof_index2(NZ)
    complex(c_double) :: A(NZA)
    complex(c_double) :: mat(NZM1*NZM2)
    integer(kint) :: conn1t(n_base1), conn2t(n_base2), i, j
    complex(kdouble) mat_t(NZM1,NZM2)

    conn1t = conn1 + 1
    conn2t = conn2 + 1

    do i = 1, NZM1
      do j = 1, NZM2
        mat_t(i,j) = mat(NZM2*(i - 1) + j)
      enddo
    enddo

    call monolis_add_matrix_to_sparse_matrix_main_C(index, item, A, n_base1, n_base2, &
      n_dof_index, n_dof_index2, conn1t, conn2t, mat_t)
  end subroutine monolis_add_matrix_to_sparse_matrix_main_C_c

  subroutine monolis_set_Dirichlet_bc_C_c(N, NZ, NZA, NZB, n_dof_index, n_dof_index2, &
    index, item, indexR, itemR, permR, A, B, node_id, ndof_bc, val) &
    bind(c, name = "monolis_set_Dirichlet_bc_C_c_main")
    implicit none
    integer(kint_c), intent(in), value :: N, NZ, NZA, NZB, node_id, ndof_bc
    integer(kint_c), intent(in) :: index(N + 1)
    integer(kint_c), intent(in) :: item(NZ)
    integer(kint_c), intent(in) :: indexR(N + 1)
    integer(kint_c), intent(in) :: itemR(NZ)
    integer(kint_c), intent(in) :: permR(NZ)
    complex(c_double), intent(in), value :: val
    integer(kint_c), intent(in) :: n_dof_index(N + 1)
    integer(kint_c), intent(in) :: n_dof_index2(NZ)
    complex(c_double) :: A(NZA)
    complex(c_double) :: B(NZB)
    integer(kint) :: nid_t, ndof_bc_t

    nid_t = node_id + 1
    ndof_bc_t = ndof_bc + 1
    call monolis_set_Dirichlet_bc_main_C(index, item, A, B, indexR, itemR, permR, &
      & n_dof_index, n_dof_index2, nid_t, ndof_bc_t, val)
  end subroutine monolis_set_Dirichlet_bc_C_c

end module mod_monolis_spmat_handler_util_wrap
