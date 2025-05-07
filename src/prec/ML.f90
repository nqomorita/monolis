!> AMG 前処理関連モジュール
module mod_monolis_precond_ML
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

  type(monolis_mat), pointer, save :: monoMAT_save
  type(monolis_com), pointer, save :: monoCOM_save

  interface
    subroutine monolis_precond_ml_setup(sym, ndof, ierr) &
        bind(c, name="monolis_ML_wrapper_setup")
      use iso_c_binding
      integer(c_int) :: sym
      integer(c_int) :: ndof
      integer(c_int) :: ierr
    end subroutine

    subroutine monolis_precond_ml_apply(rhs, ierr) &
        bind(c, name="monolis_ML_wrapper_apply")
      use iso_c_binding
      real(c_double) :: rhs(:)
      integer(c_int) :: ierr
    end subroutine

    subroutine monolis_precond_ml_clear(ierr) &
        bind(c, name="monolis_ML_wrapper_clear")
      use iso_c_binding
      integer(c_int) :: ierr
    end subroutine
  end interface

contains

  !> @ingroup prec
  !> 前処理生成：ML 前処理（実数型）
  subroutine monolis_precond_ml_setup_R(monoPRM, monoCOM, monoMAT)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in), target :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in), target :: monoMAT
    integer(kint) :: sym, ndof, ierr

    monoCOM_save => monoCOM
    monoMAT_save => monoMAT

    call monolis_std_debug_log_header("monolis_precond_ml_setup_R")
    sym = 0
    ndof = monoMAT%ndof
    call monolis_precond_ml_setup(sym, ndof, ierr);
  end subroutine monolis_precond_ml_setup_R

  !> @ingroup prec
  !> 前処理適用：SOR 前処理（実数型）
  subroutine monolis_precond_ml_apply_R(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    real(kdouble) :: X(:), Y(:)
    integer(kint) :: ierr, i

    call monolis_std_debug_log_header("monolis_precond_ml_apply_R")

    do i = 1, monoMAT%N * monoMAT%NDOF
      Y(i) = X(i)
    enddo

    call monolis_precond_ml_apply(Y, ierr)
  end subroutine monolis_precond_ml_apply_R

  !> @ingroup prec
  !> 前処理初期化：SOR 前処理（実数型）
  subroutine monolis_precond_ml_clear_R(monoPRM, monoCOM, monoMAT)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    integer(kint) :: ierr

    call monolis_std_debug_log_header("monolis_precond_ml_clear_R")

    call monolis_precond_ml_clear(ierr)
  end subroutine monolis_precond_ml_clear_R
end module mod_monolis_precond_ML

  !> wrapper section
  subroutine monolis_ML_get_nlocal(nlocal, nlocal_allcolumns, ierr) &
    bind(c, name = "monolis_ml_get_nlocal_c")
    use iso_c_binding
    use mod_monolis_utils
    use mod_monolis_precond_ML
    implicit none
    integer(c_int), intent(out) :: nlocal
    integer(c_int), intent(out) :: nlocal_allcolumns
    integer(c_int), intent(out) :: ierr
    
    nlocal = monoMAT_save%N * monoMAT_save%NDOF
    nlocal_allcolumns = monoMAT_save%NP * monoMAT_save%NDOF
    ierr = 0
  end subroutine monolis_ML_get_nlocal

  subroutine monolis_ML_getrow_nn(n_requested_rows, requested_rows, &
    allocated_space, cols, values, row_lengths, ierr) &
    bind(c, name = "monolis_ml_getrow_nn_c")
    use iso_c_binding
    use mod_monolis_utils
    use mod_monolis_precond_ML
    implicit none
    integer(c_int), intent(in) :: n_requested_rows
    integer(c_int), intent(in) :: requested_rows(n_requested_rows)
    integer(c_int), intent(in) :: allocated_space
    integer(c_int), intent(out) :: cols(allocated_space)
    real(c_double), intent(out) :: values(allocated_space)
    integer(c_int), intent(out) :: row_lengths(n_requested_rows)
    integer(c_int), intent(out) :: ierr
    integer(c_int) :: m, i, row, inod, idof, nl, nd, nu, js, je, j, jj, jdof, start, ndof
    integer(c_int) :: n

    ndof = monoMAT_save%NDOF
    m = 1
    ierr = 0
    do i = 1, n_requested_rows
      row = requested_rows(i) + 1 ! '+1' for Fortran-numbering
      inod = (row-1)/ndof + 1
      idof = row - (inod-1)*ndof
      
      n = ndof*(monoMAT_save%CSR%index(inod + 1) - monoMAT_save%CSR%index(inod))
      if (allocated_space < n) return

      start = m
      js = monoMAT_save%CSR%index(inod) + 1
      je = monoMAT_save%CSR%index(inod + 1)
      do j = js, je
        jj = monoMAT_save%CSR%item(j)
        do jdof = 1, ndof
          cols(m) = (jj-1)*ndof + jdof - 1 ! '-1' for C-numbering
          values(m) = monoMAT_save%R%A((j-1)*ndof*ndof + (idof-1)*ndof + jdof)
          m = m + 1
        enddo
      enddo
      row_lengths(i) = m - start
    enddo
    ierr = 1
  end subroutine monolis_ML_getrow_nn

  subroutine monolis_ML_matvec_nn(in_length, X, out_length, Y, ierr) &
    bind(c, name = "monolis_ml_matvec_nn_c")
    use iso_c_binding
    use mod_monolis_utils
    use mod_monolis_matvec
    use mod_monolis_precond_ML
    implicit none
    integer(c_int), intent(in) :: in_length
    real(c_double), intent(in) :: X(in_length)
    integer(c_int), intent(in) :: out_length
    real(c_double), intent(out) :: Y(out_length)
    integer(c_int), intent(out) :: ierr
    real(c_double), allocatable :: W(:)
    integer(c_int) :: i
    real(c_double) :: tspmv, tcomm

    allocate(W(monoMAT_save%NP*monoMAT_save%NDOF), source = 0.0d0)

    do i = 1, monoMAT_save%N*monoMAT_save%NDOF
      W(i) = X(i)
    enddo

    call monolis_matvec_product_main_R(monoCOM_save, monoMAT_save, W, Y, tspmv, tcomm)
    
    deallocate(W)
    ierr = 0
  end subroutine monolis_ML_matvec_nn

  subroutine monolis_ML_comm_nn(x, ierr) &
    bind(c, name = "monolis_ml_comm_nn_c")
    use iso_c_binding
    use mod_monolis_utils
    use mod_monolis_precond_ML
    implicit none
    real(c_double), intent(inout) :: x(*)
    integer(c_int), intent(out) :: ierr
    real(c_double) :: tcomm

    call monolis_mpi_update_R(monoCOM_save, monoMAT_save%ndof, x(1:monoMAT_save%NP*monoMAT_save%NDOF), tcomm)
    ierr = 0
  end subroutine monolis_ML_comm_nn
