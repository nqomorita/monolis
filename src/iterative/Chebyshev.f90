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
    integer(kint) :: i, iter, degree
    real(kdouble) :: R2, B2, eig_ratio, lambda_max
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

    eig_ratio = 30.0d0

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    !call monolis_alloc_R_1d(R, NPNDOF)
    !call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    !call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    !if(is_converge) return

    if(monoPRM%Iarray(monolis_prm_I_is_solv_prepared) == 0)then
      call monolis_palloc_R_1d(monoMAT%R%D, NPNDOF)
      call monolis_solver_JACOBI_setup(monoMAT)
      monoPRM%Iarray(monolis_prm_I_is_solv_prepared) = 1

      do i = 1, NNDOF
        monoMAT%R%D(i) = monoMAT%R%D(i) / 1.0d0
      enddo
    endif

    lambda_max = monoPRM%Rarray(monolis_prm_I_CHEBYSHEV_max_eigen_value)

    if(lambda_max == 0.0d0)then
      call monolis_solver_power_method_eigenvalue_estimation( &
        monoCOM, monoMAT, NNDOF, NPNDOF, lambda_max)
      monoPRM%Rarray(monolis_prm_I_CHEBYSHEV_max_eigen_value) = lambda_max
    endif

    !call monolis_inner_product_main_R(monoCOM, NNDOF, B, B, B2, tdotp, tcomm_dotp)

    !> Deflatin 前処理専用の最適化
    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_solver_Chebyshev_iteration( &
        monoCOM, monoMAT, NNDOF, NPNDOF, X, B, lambda_max, eig_ratio, degree, tspmv, tcomm_spmv)
      !call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      !call monolis_inner_product_main_R(monoCOM, NNDOF, R, R, R2, tdotp, tcomm_dotp)
      !call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      !if(is_converge) exit
    enddo

    call monolis_mpi_update_R_wrapper(monoCOM, monoMAT%NDOF, monoMAT%n_dof_index, X, tcomm_spmv)

    if(monoPRM%Iarray(monolis_prm_I_is_solv_prepared) == 0)then
      call monolis_dealloc_R_1d(R)
    endif
  end subroutine monolis_solver_Chebyshev

  subroutine monolis_solver_Chebyshev_iteration( &
    monoCOM, monoMAT, NNDOF, NPNDOF, X, B, lambda_max, eig_ratio, degree, tspmv, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: NNDOF, NPNDOF, degree
    real(kdouble) :: X(:), B(:)
    real(kdouble) :: lambda_max
    real(kdouble) :: tspmv, tcomm
    integer(kint) :: i, iter
    real(kdouble) :: alpha, beta, delta, theta, inv_theta
    real(kdouble) :: s1, d1, d2, rhok, rhokp1, eig_ratio
    real(kdouble), allocatable :: V(:), W(:)
    real(kdouble), pointer, contiguous :: invD(:)

    invD => monoMAT%R%D
    if(degree <= 0) return

    call monolis_alloc_R_1d(V, NPNDOF)
    call monolis_alloc_R_1d(W, NPNDOF)

    alpha = lambda_max / eig_ratio
    beta  = 1.1d0 * lambda_max
    delta = 2.0d0 / (beta - alpha)
    theta = 0.5d0 * (beta + alpha)
    s1    = theta * delta
    inv_theta = 1.0d0 / theta

    do i = 1, NNDOF
      W(i) = invD(i) * B(i) * inv_theta
      X(i) = W(i)
    enddo

    rhok = 1.0d0 / s1

    do iter = 1, degree - 1
      call monolis_matvec_product_main_R(monoCOM, monoMAT, X, V, tspmv, tcomm)
    
      rhokp1 = 1.0d0 / (2.0d0*s1 - rhok)
      d1     = rhokp1 * rhok
      d2     = 2.0d0 * rhokp1 * delta
      rhok   = rhokp1

      do i = 1, NNDOF
        W(i) = d1*W(i) + d2*invD(i)*(B(i) - V(i))
      enddo

      do i = 1, NNDOF
        X(i) = X(i) + W(i)
      enddo
    enddo

    call monolis_dealloc_R_1d(V)
    call monolis_dealloc_R_1d(W)
  end subroutine monolis_solver_Chebyshev_iteration

  subroutine monolis_solver_power_method_eigenvalue_estimation( &
    monoCOM, monoMAT, NNDOF, NPNDOF, lambda_max)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: NNDOF, NPNDOF
    real(kdouble) :: lambda_max
    real(kdouble) :: tspmv, tcomm
    integer(kint) :: iter, max_iter
    integer(kint) :: i
    real(kdouble) :: lambda, lambda_old, tol, norm_v
    real(kdouble), allocatable :: v(:), y(:)
    real(kdouble), pointer, contiguous :: invD(:)

    max_iter = 50
    tol = 1.0d-4

    invD => monoMAT%R%D

    call monolis_alloc_R_1d(v, NPNDOF)
    call monolis_alloc_R_1d(y, NPNDOF)

    call random_seed()
    call random_number(v)

    norm_v = 0.0d0
    do i = 1, NNDOF
      norm_v = norm_v + v(i)**2
    enddo

    call monolis_allreduce_R1(norm_v, monolis_mpi_sum, monoCOM%comm)
    norm_v = 1.0d0/dsqrt(norm_v)
    
    do i = 1, NNDOF
      v(i) = v(i)*norm_v
    enddo

    lambda_max = 0.0d0

    do iter = 1, max_iter
      lambda_old = lambda_max

      call monolis_matvec_product_main_R(monoCOM, monoMAT, v, y, tspmv, tcomm)

      do i = 1, NNDOF
        y(i) = invD(i) * y(i)
      enddo

      lambda = 0.0d0
      norm_v = 0.0d0
      do i = 1, NNDOF
        lambda = lambda + v(i) * y(i)
        norm_v = norm_v + v(i)**2
      enddo

      call monolis_allreduce_R1(lambda, monolis_mpi_sum, monoCOM%comm)
      call monolis_allreduce_R1(norm_v, monolis_mpi_sum, monoCOM%comm)

      lambda_max = lambda / norm_v

      norm_v = 0.0d0
      do i = 1, NNDOF
        norm_v = norm_v + y(i)**2
      enddo
      call monolis_allreduce_R1(norm_v, monolis_mpi_sum, monoCOM%comm)
      norm_v = 1.0d0/dsqrt(norm_v)

      do i = 1, NNDOF
        v(i) = y(i)*norm_v
      enddo

      if(iter > 1)then
        if(dabs(lambda_max - lambda_old) < tol * dabs(lambda_max)) exit
      endif
    enddo

    call monolis_dealloc_R_1d(v)
    call monolis_dealloc_R_1d(y)
  end subroutine monolis_solver_power_method_eigenvalue_estimation

end module mod_monolis_solver_Chebyshev
