!> LU 分解（nxn ブロック）
module mod_monolis_fact_LU_nn
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_fact_fillin
  use mod_monolis_fact_analysis
  use mod_monolis_fact_factorize
  use mod_monolis_spmat_reorder

  implicit none

  integer(kint) :: n_fact_array
  integer(kint) :: n_super_node
  integer(kint), allocatable :: super_node(:)
  integer(kint), allocatable :: fact_order(:)
  integer(kint), allocatable :: fact_array_index(:)
  integer(kint), allocatable :: add_location(:)
  real(kdouble), allocatable :: fact_array(:)

contains

  !> @ingroup prec
  !> 前処理生成：LU 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_LU_nn_setup_R(monoMAT, monoPREC)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    type(monolis_mat) :: monoMAT_reorder
    logical :: is_asym = .false.
    logical :: is_fillin = .true.
    real(kdouble) :: t(100)

    t(1) = monolis_get_time_global_sync()

    !> analysis phase
    !call monolis_matrix_reordering_fw_R(monoMAT, monoMAT_reorder)

    t(2) = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"monolis_matrix_reordering_fw_R             ", t(2) - t(1)

    call monolis_matrix_get_fillin(monoMAT, monoPREC, is_asym, is_fillin)
    !call monolis_matrix_get_fillin(monoMAT_reorder, monoPREC, is_asym, is_fillin)

    t(3) = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"monolis_matrix_get_fillin                  ", t(3) - t(2)

    call monolis_matrix_alloc_with_fillin(monoPREC, is_asym)

    t(4) = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"monolis_matrix_alloc_with_fillin           ", t(4) - t(3)

    call monolis_matrix_get_super_node_information(monoPREC, n_super_node, super_node)

    t(5) = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"monolis_matrix_get_super_node_information  ", t(5) - t(4)

    call monolis_matrix_get_factorize_order(monoPREC, fact_order)

    t(6) = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"monolis_matrix_get_factorize_order         ", t(6) - t(5)

    call monolis_matrix_get_factorize_array(monoPREC, fact_order, n_fact_array, fact_array, fact_array_index)

    t(7) = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"monolis_matrix_get_factorize_array         ", t(7) - t(6)

    call monolis_matrix_set_value_of_factorize_array(monoMAT, monoPREC, &
    !call monolis_matrix_set_value_of_factorize_array(monoMAT_reorder, monoPREC, &
      & fact_order, n_fact_array, fact_array, fact_array_index)

    t(8) = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"monolis_matrix_set_value_of_factorize_array", t(8) - t(7)

    call monolis_matrix_get_add_location(monoPREC, fact_order, fact_array_index, add_location)

    t(9) = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"monolis_matrix_get_add_location            ", t(9) - t(8)

    !> factorization phase
    call monolis_matrix_factorize_mf(monoPREC, fact_order, fact_array, fact_array_index, add_location)

    t(10) = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"monolis_matrix_factorize_mf                ", t(10) - t(9)

    call monolis_matrix_copy_lu_factor(monoPREC, fact_order, fact_array, fact_array_index)

    t(11) = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"monolis_matrix_copy_lu_factor", t(11) - t(10)
  end subroutine monolis_fact_LU_nn_setup_R

  !> @ingroup prec
  !> 前処理生成：LU 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_LU_nn_setup_C(monoMAT, monoPREC)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC

  end subroutine monolis_fact_LU_nn_setup_C

  !> @ingroup prec
  !> 前処理適用：LU 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_LU_nn_apply_R(monoMAT, monoPREC, Y, X)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in] 前処理構造体
    type(monolis_mat), target, intent(in) :: monoPREC
    real(kdouble) :: X(:), Y(:)
    integer(kint) :: N, NDOF
    integer(kint) :: i, j, k, in, jS, jE, kn
    real(kdouble) :: X1, A1
    real(kdouble), allocatable :: S(:)
    integer(kint), pointer :: idxU(:), itemU(:)
    real(kdouble), pointer :: A(:)

    N  = monoPREC%N
    NDOF = 1
    allocate(S(N), source = 0.0d0)

    !S = Y(1:N)
    !call monolis_reorder_vector_fw(monoMAT, N, NDOF, S, X)
    X = Y(1:N)

    idxU => monoPREC%SCSR%indexU
    itemU => monoPREC%SCSR%itemU
    A => monoPREC%R%A
    S = 0.0d0

    !L
    do i = 1, N
      A1 = S(i)
      in = idxU(i) + 1
      X(i) = A(in)*(X(i) - A1)
      X1 = X(i)
      jS = idxU(i) + 2
      jE = idxU(i + 1)
      do j = jS, jE
        in = itemU(j)
        S(in) = S(in) + A(j)*X1
      enddo
    enddo

    !D
    do i = 1, N
      in = idxU(i) + 1
      X(i) = X(i)/A(in)
    enddo

    !U
    do i = N, 1, -1
      A1 = 0.0d0
      jS = idxU(i) + 2
      jE = idxU(i + 1)
      do j = jE, jS, -1
        in = itemU(j)
        X1 = X(in)
        A1 = A1 + A(j)*X1
      enddo
      in = idxU(i) + 1
      X(i) = A(in)*(X(i) - A1)
    enddo

    !S = X(1:N)
    !call monolis_reorder_back_vector_bk(monoMAT, N, NDOF, S, X)
    deallocate(S)
  end subroutine monolis_fact_LU_nn_apply_R

  !> 前処理適用：LU 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_LU_nn_apply_C(monoMAT, monoPREC, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    complex(kdouble) :: X(:), Y(:)

  end subroutine monolis_fact_LU_nn_apply_C

  !> @ingroup prec
  !> 前処理初期化：LU 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_LU_nn_clear_R(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_R_1d(monoPREC%R%D)
  end subroutine monolis_fact_LU_nn_clear_R

  !> @ingroup prec
  !> 前処理初期化：LU 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_LU_nn_clear_C(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_C_1d(monoPREC%C%D)
  end subroutine monolis_fact_LU_nn_clear_C
end module mod_monolis_fact_LU_nn
