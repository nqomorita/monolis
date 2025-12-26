!> Chebyshev polynomial smoother module
module mod_monolis_solver_Chebyshev
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_inner_product
  use mod_monolis_vec_util
  use mod_monolis_converge
  use mod_monolis_solver_JACOBI

  implicit none

  private
  public monolis_solver_Chebyshev

contains

  subroutine monolis_solver_Chebyshev(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] parameter structure
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] communication table structure
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] matrix structure
    type(monolis_mat), target, intent(inout) :: monoMAT
    !> [in,out] preconditioner structure
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: NNDOF, NPNDOF
    integer(kint) :: iter, degree
    real(kdouble) :: R2, B2, lambda_min, lambda_max
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    real(kdouble), pointer, contiguous :: B(:), X(:)
    real(kdouble), allocatable :: R(:)
    logical :: is_converge

    X => monoMAT%R%X
    B => monoMAT%R%B

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    degree = monoPRM%Iarray(monolis_prm_I_CHEBYSHEV_degree)
    if(degree <= 0) degree = 3

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    !call monolis_alloc_R_1d(R, NPNDOF)
    !call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    !call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    !if(is_converge) return

    lambda_min = monoPRM%Rarray(monolis_prm_I_CHEBYSHEV_min_eigen_value)
    lambda_max = monoPRM%Rarray(monolis_prm_I_CHEBYSHEV_max_eigen_value)

    if(lambda_min == 0.0d0 .or. lambda_max == 0.0d0)then
      call monolis_solver_power_method_eigenvalue_estimation( &
        monoCOM, monoMAT, NNDOF, NPNDOF, lambda_min, lambda_max)
      monoPRM%Rarray(monolis_prm_I_CHEBYSHEV_min_eigen_value) = lambda_min
      monoPRM%Rarray(monolis_prm_I_CHEBYSHEV_max_eigen_value) = lambda_max
    endif

    !call monolis_inner_product_main_R(monoCOM, NNDOF, B, B, B2, tdotp, tcomm_dotp)

    !> Deflatin 前処理専用の最適化
    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_solver_Chebyshev_iteration( &
        monoCOM, monoMAT, NNDOF, NPNDOF, X, B, lambda_min, lambda_max, degree, tspmv, tcomm_spmv)
      !call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      !call monolis_inner_product_main_R(monoCOM, NNDOF, R, R, R2, tdotp, tcomm_dotp)
      !call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      !if(is_converge) exit
    enddo

    call monolis_mpi_update_R_wrapper(monoCOM, monoMAT%NDOF, monoMAT%n_dof_index, X, tcomm_spmv)

    call monolis_dealloc_R_1d(R)
  end subroutine monolis_solver_Chebyshev

  subroutine monolis_solver_Chebyshev_iteration( &
    monoCOM, monoMAT, NNDOF, NPNDOF, X, B, lambda_min, lambda_max, degree, tspmv, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: NNDOF, NPNDOF, degree
    real(kdouble) :: X(:), B(:)
    real(kdouble) :: lambda_min, lambda_max
    real(kdouble) :: tspmv, tcomm
    integer(kint) :: k, i
    real(kdouble) :: alpha, beta, c, d, theta
    real(kdouble) :: alpha_new, beta_new
    real(kdouble), allocatable :: P(:), R(:)

    call monolis_alloc_R_1d(P, NPNDOF)
    call monolis_alloc_R_1d(R, NPNDOF)

    ! Initialize Chebyshev polynomial parameters
    c = (lambda_max - lambda_min) / 2.0d0
    d = (lambda_max + lambda_min) / 2.0d0
    alpha = 0.0d0

    do k = 1, degree
      call monolis_matvec_product_main_R(monoCOM, monoMAT, X, R, tspmv, tcomm)
      do i = 1, NNDOF
        R(i) = B(i) - R(i)
      enddo

      if(k == 1)then
        alpha = 1.0d0 / d
      elseif(k == 2)then
        alpha = 2.0d0*d / (2.0d0*d*d - c*c)
      elseif(k >= 3)then
        alpha = 4.0d0 / (d - alpha*c*c)
      endif

      beta = alpha*d - 1.0d0

      do i = 1, NNDOF
        P(i) = alpha*R(i) - beta*P(i)
        X(i) = X(i) + P(i)
      enddo
    enddo

    call monolis_dealloc_R_1d(P)
    call monolis_dealloc_R_1d(R)
  end subroutine monolis_solver_Chebyshev_iteration

end module mod_monolis_solver_Chebyshev
