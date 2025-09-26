!> Chebyshev polynomial smoother module
module mod_monolis_solver_CHEBYSHEV
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_inner_product
  use mod_monolis_vec_util
  use mod_monolis_converge

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
    real(kdouble), pointer :: B(:), X(:)
    real(kdouble), allocatable :: R(:)
    logical :: is_converge

    X => monoMAT%R%X
    B => monoMAT%R%B

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    ! Get Chebyshev polynomial degree (default: 3)
    degree = monoPRM%Iarray(monolis_prm_I_precond_degree)
    if(degree <= 0) degree = 3

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    call monolis_alloc_R_1d(R, NPNDOF)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)

    if(is_converge) return

    call monolis_solver_CHEBYSHEV_eigenvalue_estimation( &
      monoCOM, monoMAT, NNDOF, NPNDOF, lambda_min, lambda_max, tspmv, tcomm_spmv)

    call monolis_inner_product_main_R(monoCOM, NNDOF, B, B, B2, tdotp, tcomm_dotp)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_solver_CHEBYSHEV_iteration( &
        monoCOM, monoMAT, NNDOF, NPNDOF, X, B, lambda_min, lambda_max, degree, tspmv, tcomm_spmv)
      call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      call monolis_inner_product_main_R(monoCOM, NNDOF, R, R, R2, tdotp, tcomm_dotp)
      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit
    enddo

    call monolis_mpi_update_R_wrapper(monoCOM, monoMAT%NDOF, monoMAT%n_dof_index, X, tcomm_spmv)

    call monolis_dealloc_R_1d(R)
  end subroutine monolis_solver_Chebyshev

  subroutine monolis_solver_CHEBYSHEV_iteration( &
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
    real(kdouble), allocatable :: P0(:), P1(:), P2(:), AX(:), R(:)

    call monolis_alloc_R_1d(P0, NPNDOF)
    call monolis_alloc_R_1d(P1, NPNDOF)
    call monolis_alloc_R_1d(P2, NPNDOF)
    call monolis_alloc_R_1d(AX, NPNDOF)
    call monolis_alloc_R_1d(R, NPNDOF)

    ! Compute initial residual: R = B - A*X
    call monolis_matvec_product_main_R(monoCOM, monoMAT, X, AX, tspmv, tcomm)
    do i = 1, NNDOF
      R(i) = B(i) - AX(i)
    enddo

    ! Initialize Chebyshev polynomial parameters
    c = (lambda_max + lambda_min) / 2.0d0
    d = (lambda_max - lambda_min) / 2.0d0
    
    ! First iteration (k=0): P0 = (2/d) * R
    alpha = 2.0d0 / d
    do i = 1, NNDOF
      P0(i) = alpha * R(i)
      X(i) = X(i) + P0(i)
    enddo

    if(degree == 1) goto 100

    ! Second iteration (k=1)
    call monolis_matvec_product_main_R(monoCOM, monoMAT, X, AX, tspmv, tcomm)
    do i = 1, NNDOF
      R(i) = B(i) - AX(i)
    enddo

    theta = d**2 / (2.0d0*c**2 - d**2)
    alpha_new = 4.0d0 * theta / d
    beta_new = 2.0d0 * theta

    do i = 1, NNDOF
      P1(i) = alpha_new * R(i) - beta_new * P0(i)
      X(i) = X(i) + P1(i)
    enddo

    if(degree == 2) goto 100

    ! Higher order iterations (k >= 2)
    do k = 2, degree - 1
      call monolis_matvec_product_main_R(monoCOM, monoMAT, X, AX, tspmv, tcomm)
      do i = 1, NNDOF
        R(i) = B(i) - AX(i)
      enddo

      theta = 4.0d0*c**2 - d**2
      alpha = 8.0d0 * c / d
      beta = 2.0d0

      do i = 1, NNDOF
        P2(i) = alpha * R(i) - beta * P1(i)
        X(i) = X(i) + P2(i)
      enddo

      ! Shift polynomials for next iteration
      P0(:) = P1(:)
      P1(:) = P2(:)
    enddo

100 continue

    call monolis_dealloc_R_1d(P0)
    call monolis_dealloc_R_1d(P1)
    call monolis_dealloc_R_1d(P2)
    call monolis_dealloc_R_1d(AX)
    call monolis_dealloc_R_1d(R)
  end subroutine monolis_solver_CHEBYSHEV_iteration

  subroutine monolis_solver_CHEBYSHEV_eigenvalue_estimation( &
    monoCOM, monoMAT, NNDOF, NPNDOF, lambda_min, lambda_max, tspmv, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: NNDOF, NPNDOF
    real(kdouble) :: lambda_min, lambda_max
    real(kdouble) :: tspmv, tcomm
    integer(kint) :: iter, max_iter
    integer(kint) :: i
    real(kdouble) :: lambda, lambda_old, tol, norm_v
    real(kdouble), allocatable :: v(:), Av(:)

    max_iter = 50
    tol = 1.0d-6

    call monolis_alloc_R_1d(v, NPNDOF)
    call monolis_alloc_R_1d(Av, NPNDOF)

    ! Initialize random vector for power method
    call random_seed()

    do i = 1, NNDOF
      call random_number(v(i))
      v(i) = v(i) - 0.5d0  ! center around zero
    enddo

    ! Normalize initial vector
    norm_v = 0.0d0
    do i = 1, NNDOF
      norm_v = norm_v + v(i)**2
    enddo
    call monolis_allreduce_R1(norm_v, monolis_mpi_sum, monoCOM%comm)
    norm_v = dsqrt(norm_v)
    
    if(norm_v > 0.0d0) then
      do i = 1, NNDOF
        v(i) = v(i) / norm_v
      enddo
    endif

    lambda_max = 0.0d0

    ! Power method for maximum eigenvalue
    do iter = 1, max_iter
      lambda_old = lambda_max

      ! Av = A * v
      call monolis_matvec_product_main_R(monoCOM, monoMAT, v, Av, tspmv, tcomm)

      ! lambda = v^T * Av / (v^T * v)
      lambda = 0.0d0
      norm_v = 0.0d0
      do i = 1, NNDOF
        lambda = lambda + v(i) * Av(i)
        norm_v = norm_v + v(i)**2
      enddo

      call monolis_allreduce_R1(lambda, monolis_mpi_sum, monoCOM%comm)
      call monolis_allreduce_R1(norm_v, monolis_mpi_sum, monoCOM%comm)

      lambda_max = lambda / norm_v

      ! Normalize Av for next iteration
      norm_v = 0.0d0
      do i = 1, NNDOF
        norm_v = norm_v + Av(i)**2
      enddo
      call monolis_allreduce_R1(norm_v, monolis_mpi_sum, monoCOM%comm)
      norm_v = dsqrt(norm_v)

      if(norm_v > 0.0d0) then
        do i = 1, NNDOF
          v(i) = Av(i) / norm_v
        enddo
      endif

      ! Check convergence
      if(iter > 1 .and. dabs(lambda_max - lambda_old) < tol * dabs(lambda_max)) exit
    enddo

    ! Estimate minimum eigenvalue using Gershgorin circle theorem
    ! For symmetric positive definite matrices, use a simple approximation
    lambda_min = 0.1d0 * lambda_max

    ! Ensure proper bounds
    if(lambda_min <= 0.0d0) lambda_min = 0.1d0
    if(lambda_max <= lambda_min) lambda_max = 10.0d0 * lambda_min

    call monolis_dealloc_R_1d(v)
    call monolis_dealloc_R_1d(Av)
  end subroutine monolis_solver_CHEBYSHEV_eigenvalue_estimation

end module mod_monolis_solver_CHEBYSHEV
