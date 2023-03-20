!> 疎行列操作関数群（メイン関数）
module mod_monolis_spmat_handler_util_wrap
  use mod_monolis_utils
  use mod_monolis_spmat_handler
  use mod_monolis_spmat_handler_util
  use iso_c_binding

  implicit none

contains

  subroutine monolis_set_scalar_to_sparse_matrix_R_c(N, NZ, NDOF, index, item, A, i, j, sub_i, sub_j, val) &
    & bind(c, name = "monolis_set_scalar_to_sparse_matrix_R_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, i, j, sub_i, sub_j
    integer(c_int), intent(in) :: index(N + 1)
    integer(c_int), intent(in) :: item(NZ)
    real(c_double) :: A(NDOF*NDOF*NZ)
    real(c_double), value :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_set_scalar_to_sparse_matrix_main_R(index, item, A, NDOF, i_t, j_t, sub_i_t, sub_j_t, val)
  end subroutine monolis_set_scalar_to_sparse_matrix_R_c

  subroutine monolis_add_scalar_to_sparse_matrix_R_c(N, NZ, NDOF, index, item, A, i, j, sub_i, sub_j, val) &
    & bind(c, name = "monolis_add_scalar_to_sparse_matrix_R_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, i, j, sub_i, sub_j
    integer(c_int), intent(in) :: index(N + 1)
    integer(c_int), intent(in) :: item(NZ)
    real(c_double) :: A(NDOF*NDOF*NZ)
    real(c_double), value :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_add_scalar_to_sparse_matrix_main_R(index, item, A, NDOF, i_t, j_t, sub_i_t, sub_j_t, val)
  end subroutine monolis_add_scalar_to_sparse_matrix_R_c

  subroutine monolis_get_scalar_from_sparse_matrix_R_c(N, NZ, NDOF, index, item, A, i, j, sub_i, sub_j, val, is_find) &
    & bind(c, name = "monolis_get_scalar_from_sparse_matrix_R_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, i, j, sub_i, sub_j
    integer(c_int), intent(in) :: index(N + 1)
    integer(c_int), intent(in) :: item(NZ)
    integer(c_int) :: is_find
    real(c_double) :: A(NDOF*NDOF*NZ)
    real(c_double) :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t
    logical :: is_find_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_get_scalar_from_sparse_matrix_main_R(index, item, A, NDOF, i_t, j_t, sub_i_t, sub_j_t, val, is_find_t)
    is_find = 0
    if(is_find_t) is_find = 1
  end subroutine monolis_get_scalar_from_sparse_matrix_R_c

  subroutine monolis_add_matrix_to_sparse_matrix_main_R_c(N, NZ, NDOF, index, item, A, &
    & n_base1, n_base2, conn1, conn2, mat) &
    & bind(c, name = "monolis_add_matrix_to_sparse_matrix_main_R_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, n_base1, n_base2
    integer(c_int), intent(in) :: index(N + 1)
    integer(c_int), intent(in) :: item(NZ)
    integer(c_int), intent(in) :: conn1(n_base1)
    integer(c_int), intent(in) :: conn2(n_base2)
    real(c_double) :: A(NDOF*NDOF*NZ)
    real(c_double) :: mat(NDOF*NDOF*n_base1*n_base2)
    integer(kint) :: conn1t(n_base1), conn2t(n_base2), i, j
    real(kdouble) mat_t(NDOF*n_base1,NDOF*n_base2)

    conn1t = conn1 + 1
    conn2t = conn2 + 1

    do i = 1, NDOF*n_base1
      do j = 1, NDOF*n_base2
        mat_t(i,j) = mat(NDOF*n_base2*(i - 1) + j)
      enddo
    enddo

    call monolis_add_matrix_to_sparse_matrix_main_R(index, item, A, n_base1, n_base2, NDOF, conn1t, conn2t, mat_t)
  end subroutine monolis_add_matrix_to_sparse_matrix_main_R_c

  subroutine monolis_set_Dirichlet_bc_R_c(N, NZ, NDOF, index, item, indexR, itemR, permR, &
    & A, B, node_id, ndof_bc, val) &
    & bind(c, name = "monolis_set_Dirichlet_bc_R_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, node_id, ndof_bc
    integer(c_int), intent(in) :: index(N + 1)
    integer(c_int), intent(in) :: item(NZ)
    integer(c_int), intent(in) :: indexR(N + 1)
    integer(c_int), intent(in) :: itemR(NZ)
    integer(c_int), intent(in) :: permR(NZ)
    real(c_double), intent(in), value :: val
    real(c_double) :: A(NDOF*NDOF*NZ)
    real(c_double) :: B(NDOF*N)
    integer(kint) :: nid_t, ndof_bc_t

    nid_t = node_id + 1
    ndof_bc_t = ndof_bc + 1
    call monolis_set_Dirichlet_bc_main_R(index, item, A, B, indexR, itemR, permR, &
      & ndof, nid_t, ndof_bc_t, val)
  end subroutine monolis_set_Dirichlet_bc_R_c

  subroutine monolis_set_scalar_to_sparse_matrix_C_c(N, NZ, NDOF, index, item, A, i, j, sub_i, sub_j, val) &
    & bind(c, name = "monolis_set_scalar_to_sparse_matrix_C_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, i, j, sub_i, sub_j
    integer(c_int), intent(in) :: index(N + 1)
    integer(c_int), intent(in) :: item(NZ)
    complex(c_double) :: A(NDOF*NDOF*NZ)
    complex(c_double), value :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_set_scalar_to_sparse_matrix_main_C(index, item, A, NDOF, i_t, j_t, sub_i_t, sub_j_t, val)
  end subroutine monolis_set_scalar_to_sparse_matrix_C_c

  subroutine monolis_add_scalar_to_sparse_matrix_C_c(N, NZ, NDOF, index, item, A, i, j, sub_i, sub_j, val) &
    & bind(c, name = "monolis_add_scalar_to_sparse_matrix_C_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, i, j, sub_i, sub_j
    integer(c_int), intent(in) :: index(N + 1)
    integer(c_int), intent(in) :: item(NZ)
    complex(c_double) :: A(NDOF*NDOF*NZ)
    complex(c_double), value :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_add_scalar_to_sparse_matrix_main_C(index, item, A, NDOF, i_t, j_t, sub_i_t, sub_j_t, val)
  end subroutine monolis_add_scalar_to_sparse_matrix_C_c

  subroutine monolis_get_scalar_from_sparse_matrix_C_c(N, NZ, NDOF, index, item, A, i, j, sub_i, sub_j, val, is_find) &
    & bind(c, name = "monolis_get_scalar_from_sparse_matrix_C_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, i, j, sub_i, sub_j
    integer(c_int), intent(in) :: index(N + 1)
    integer(c_int), intent(in) :: item(NZ)
    integer(c_int) :: is_find
    complex(c_double) :: A(NDOF*NDOF*NZ)
    complex(c_double) :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t
    logical :: is_find_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_get_scalar_from_sparse_matrix_main_C(index, item, A, NDOF, i_t, j_t, sub_i_t, sub_j_t, val, is_find_t)
    is_find = 0
    if(is_find_t) is_find = 1
  end subroutine monolis_get_scalar_from_sparse_matrix_C_c

  subroutine monolis_add_matrix_to_sparse_matrix_main_C_c(N, NZ, NDOF, index, item, A, &
    & n_base1, n_base2, conn1, conn2, mat) &
    & bind(c, name = "monolis_add_matrix_to_sparse_matrix_main_C_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, n_base1, n_base2
    integer(c_int), intent(in) :: index(N + 1)
    integer(c_int), intent(in) :: item(NZ)
    integer(c_int), intent(in) :: conn1(n_base1)
    integer(c_int), intent(in) :: conn2(n_base2)
    complex(c_double) :: A(NDOF*NDOF*NZ)
    complex(c_double) :: mat(NDOF*NDOF*n_base1*n_base2)
    integer(kint) :: conn1t(n_base1), conn2t(n_base2), i, j
    complex(kdouble) mat_t(NDOF*n_base1,NDOF*n_base2)

    conn1t = conn1 + 1
    conn2t = conn2 + 1

    do i = 1, NDOF*n_base1
      do j = 1, NDOF*n_base2
        mat_t(i,j) = mat(NDOF*n_base2*(i - 1) + j)
      enddo
    enddo

    call monolis_add_matrix_to_sparse_matrix_main_C(index, item, A, n_base1, n_base2, NDOF, conn1t, conn2t, mat_t)
  end subroutine monolis_add_matrix_to_sparse_matrix_main_C_c

  subroutine monolis_set_Dirichlet_bc_C_c(N, NZ, NDOF, index, item, indexR, itemR, permR, &
    & A, B, node_id, ndof_bc, val) &
    & bind(c, name = "monolis_set_Dirichlet_bc_C_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, node_id, ndof_bc
    integer(c_int), intent(in) :: index(N + 1)
    integer(c_int), intent(in) :: item(NZ)
    integer(c_int), intent(in) :: indexR(N + 1)
    integer(c_int), intent(in) :: itemR(NZ)
    integer(c_int), intent(in) :: permR(NZ)
    complex(c_double), intent(in), value :: val
    complex(c_double) :: A(NDOF*NDOF*NZ)
    complex(c_double) :: B(NDOF*N)
    integer(kint) :: nid_t, ndof_bc_t

    nid_t = node_id + 1
    ndof_bc_t = ndof_bc + 1
    call monolis_set_Dirichlet_bc_main_C(index, item, A, B, indexR, itemR, permR, &
      & ndof, nid_t, ndof_bc_t, val)
  end subroutine monolis_set_Dirichlet_bc_C_c

end module mod_monolis_spmat_handler_util_wrap
