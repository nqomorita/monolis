!> 固有値ソルバテストモジュール
module mod_monolis_eigen_solver_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_eigen_solver_test()
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: mat
    !> monolis 通信構造体
    type(monolis_com) :: com
    !> 取得固有値数
    integer(kint) :: n_get_eigen
    !> 収束判定閾値
    real(kdouble) :: ths
    !> 最大反復回数
    integer(kint) :: maxiter
    !> 固有値
    real(kdouble) :: eig_val(6)
    !> 固有ベクトル
    real(kdouble) :: eig_mode(6,6)
    !> Dirhchlet 境界条件判定フラグ
    logical :: is_bc(6)
    integer(kint) :: nnode, nelem, elem(2,9), n_dof
    integer(kint) :: i
    real(kdouble) :: r_ans(6)

    call monolis_std_global_log_string("lanczos_initialze")
    call monolis_std_global_log_string("monolis_gram_schmidt_R")
    call monolis_std_global_log_string("monolis_eigen_standard_lanczos_R")
    call monolis_std_global_log_string("monolis_eigen_standard_lanczos_R_main")
    call monolis_std_global_log_string("monolis_eigen_inverted_standard_lanczos_R")
    call monolis_std_global_log_string("monolis_eigen_inverted_standard_lanczos_R_main")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    nnode = 6

    nelem = 5

    n_dof = 1

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;
    elem(1,4) = 4; elem(2,4) = 5;
    elem(1,5) = 5; elem(2,5) = 6;

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, nnode, 2, n_dof, nelem, elem)

    do i = 1, 6
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, 1, 1, 2.0d0)
    enddo

    do i = 1, 4
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i + 1, 1, 1, 1.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i + 1, i, 1, 1, 1.0d0)
    enddo

    n_get_eigen = 5

    is_bc = .false.
    is_bc(6) = .true.

    ths = 1.0d0-8

    maxiter = 5

    eig_val = 0.0d0

    eig_mode = 0.0d0

    call monolis_eigen_standard_lanczos_R &
      & (mat, com, n_get_eigen, ths, maxiter, eig_val, eig_mode, is_bc)

    call monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 1 a", eig_val(5), 0.267949192431122d0)
    call monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 2 a", eig_val(4), 1.0d0)
    call monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 3 a", eig_val(3), 2.0d0)
    call monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 4 a", eig_val(2), 3.0d0)
    call monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 5 a", eig_val(1), 3.732050807568877d0)

    r_ans(1) =-0.28867513459481281d0
    r_ans(2) = 0.5d0
    r_ans(3) =-0.57735026918962640d0
    r_ans(4) = 0.5d0
    r_ans(5) =-0.28867513459481270d0
    r_ans(6) = 0.0d0

    call monolis_test_check_eq_R("monolis_eigen_standard_lanczos_R 1 b", dabs(r_ans), dabs(eig_mode(:,5)))

    r_ans(1) = 0.5d0
    r_ans(2) =-0.5d0
    r_ans(3) = 0.0d0
    r_ans(4) = 0.5d0
    r_ans(5) =-0.5d0
    r_ans(6) = 0.0d0

    call monolis_test_check_eq_R("monolis_eigen_standard_lanczos_R 2 b", dabs(r_ans), dabs(eig_mode(:,4)))

    r_ans(1) =-0.5773502691896257d0
    r_ans(2) = 0.0d0
    r_ans(3) = 0.5773502691896257d0
    r_ans(4) = 0.0d0
    r_ans(5) =-0.5773502691896257d0
    r_ans(6) = 0.0d0

    call monolis_test_check_eq_R("monolis_eigen_standard_lanczos_R 3 b", dabs(r_ans), dabs(eig_mode(:,3)))

    r_ans(1) = 0.5d0
    r_ans(2) = 0.5d0
    r_ans(3) = 0.0d0
    r_ans(4) =-0.5d0
    r_ans(5) =-0.5d0
    r_ans(6) = 0.0d0

    call monolis_test_check_eq_R("monolis_eigen_standard_lanczos_R 4 b", dabs(r_ans), dabs(eig_mode(:,2)))

    r_ans(1) =-0.28867513459481281d0
    r_ans(2) =-0.5d0
    r_ans(3) =-0.57735026918962640d0
    r_ans(4) =-0.5d0
    r_ans(5) =-0.28867513459481270d0
    r_ans(6) = 0.0d0

    call monolis_test_check_eq_R("monolis_eigen_standard_lanczos_R 5 b", dabs(r_ans), dabs(eig_mode(:,1)))

    eig_val = 0.0d0

    eig_mode = 0.0d0

    call monolis_eigen_inverted_standard_lanczos_R &
      & (mat, com, n_get_eigen, ths, maxiter, eig_val, eig_mode, is_bc)

    call monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 1 c", eig_val(1), 0.267949192431122d0)
    call monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 2 c", eig_val(2), 1.0d0)
    call monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 3 c", eig_val(3), 2.0d0)
    call monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 4 c", eig_val(4), 3.0d0)
    call monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 5 c", eig_val(5), 3.732050807568877d0)

    r_ans(1) = 0.28867513459481281d0
    r_ans(2) =-0.5d0
    r_ans(3) = 0.57735026918962640d0
    r_ans(4) =-0.5d0
    r_ans(5) = 0.28867513459481270d0
    r_ans(6) = 0.0d0

    call monolis_test_check_eq_R("monolis_eigen_standard_lanczos_R 1 d", dabs(r_ans), dabs(eig_mode(:,1)))

    r_ans(1) = 0.5d0
    r_ans(2) =-0.5d0
    r_ans(3) = 0.0d0
    r_ans(4) = 0.5d0
    r_ans(5) =-0.5d0
    r_ans(6) = 0.0d0

    call monolis_test_check_eq_R("monolis_eigen_standard_lanczos_R 2 d", dabs(r_ans), dabs(eig_mode(:,2)))

    r_ans(1) = 0.5773502691896257d0
    r_ans(2) = 0.0d0
    r_ans(3) =-0.5773502691896257d0
    r_ans(4) = 0.0d0
    r_ans(5) = 0.5773502691896257d0
    r_ans(6) = 0.0d0

    call monolis_test_check_eq_R("monolis_eigen_standard_lanczos_R 3 d", dabs(r_ans), dabs(eig_mode(:,3)))

    r_ans(1) = 0.5d0
    r_ans(2) = 0.5d0
    r_ans(3) = 0.0d0
    r_ans(4) =-0.5d0
    r_ans(5) =-0.5d0
    r_ans(6) = 0.0d0

    call monolis_test_check_eq_R("monolis_eigen_standard_lanczos_R 4 d", dabs(r_ans), dabs(eig_mode(:,4)))

    r_ans(1) =-0.28867513459481281d0
    r_ans(2) =-0.5d0
    r_ans(3) =-0.57735026918962640d0
    r_ans(4) =-0.5d0
    r_ans(5) =-0.28867513459481270d0
    r_ans(6) = 0.0d0

    call monolis_test_check_eq_R("monolis_eigen_standard_lanczos_R 5 d", dabs(r_ans), dabs(eig_mode(:,5)))

  end subroutine monolis_eigen_solver_test
end module mod_monolis_eigen_solver_test

