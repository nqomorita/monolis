!> JACOBI 法モジュール
module mod_monolis_solver_JACOBI
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
  public monolis_solver_JACOBI
  public monolis_solver_power_method_eigenvalue_estimation

contains

  subroutine monolis_solver_JACOBI(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), target, intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: NNDOF, NPNDOF
    integer(kint) :: iter
    real(kdouble) :: R2, B2, omega, no_use
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    real(kdouble), pointer, contiguous :: B(:), X(:)
    !real(kdouble), allocatable :: R(:)
    logical :: is_converge

    X => monoMAT%R%X
    B => monoMAT%R%B

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)
    omega = monoPRM%Rarray(monolis_prm_R_DCG_inner_relaxation_factor)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

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
    endif
    !call monolis_inner_product_main_R(monoCOM, NNDOF, B, B, B2, tdotp, tcomm_dotp)

    if(omega <= 0.0d0)then
      call monolis_solver_JACOBI_get_auto_relax_factor(monoCOM, monoMAT, NPNDOF, omega)
      monoPRM%Rarray(monolis_prm_R_DCG_inner_relaxation_factor) = omega
      !call monolis_solver_power_method_eigenvalue_estimation( &
      !  monoCOM, monoMAT, NNDOF, NPNDOF, no_use, omega)
      !monoPRM%Rarray(monolis_prm_R_DCG_inner_relaxation_factor) = omega
    endif

    !> Deflatin 前処理専用の最適化
    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_solver_JACOBI_matvec(monoCOM, monoMAT, NNDOF, NPNDOF, X, B, omega, tspmv, tcomm_spmv)
      !call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      !call monolis_inner_product_main_R(monoCOM, NNDOF, R, R, R2, tdotp, tcomm_dotp)
      !call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      !if(is_converge) exit
    enddo

    call monolis_mpi_update_R_wrapper(monoCOM, monoMAT%NDOF, monoMAT%n_dof_index, X, tcomm_spmv)

    !call monolis_dealloc_R_1d(R)
    if(monoPRM%Iarray(monolis_prm_I_is_solv_prepared) == 0)then
      call monolis_pdealloc_R_1d(monoMAT%R%D)
    endif
  end subroutine monolis_solver_JACOBI

  subroutine monolis_solver_JACOBI_setup(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, ii, j, jS, jE, in, N, n1, n2, nz
    integer(kint), pointer :: index(:), item(:)
    real(kdouble), pointer :: A(:), D(:)

    A => monoMAT%R%A
    D => monoMAT%R%D
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    do i = 1, monoMAT%N
      jS = index(i) + 1
      jE = index(i + 1)
      n1 = monoMAT%n_dof_list(i)
      n2 = monoMAT%n_dof_index(i)
      do ii = jS, jE
        in = item(ii)
        nz = monoMAT%n_dof_index2(ii)
        if(i == in)then
          do j = 1, n1
            D(n2 + j) = A(nz + (n1+1)*(j-1) + 1)
          enddo
        endif
      enddo
    enddo
  end subroutine monolis_solver_JACOBI_setup

  subroutine monolis_solver_JACOBI_get_auto_relax_factor(monoCOM, monoMAT, NPNDOF, omega)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat), target :: monoMAT
    integer(kint) :: NPNDOF
    integer(kint) :: i, ii, j, jS, jE, in, n1, n2, n3, nz, dof_i, dof_j
    real(kdouble) :: omega, U, ri, ri_max, s_i, s_j, Aij_abs
    real(kdouble) :: tcomm_spmv
    real(kdouble), allocatable :: s(:)
    integer(kint), pointer :: index(:), item(:)
    real(kdouble), pointer :: A(:), D(:)

    A => monoMAT%R%A
    D => monoMAT%R%D
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    call monolis_mpi_update_R_wrapper(monoCOM, monoMAT%NDOF, monoMAT%n_dof_index, D, tcomm_spmv)

    !# s[i] = 1/sqrt(Aii) を計算
    call monolis_alloc_R_1d(s, NPNDOF)
    
    do i = 1, NPNDOF
      s(i) = 1.0d0 / dsqrt(dabs(D(i)))
    enddo

    !# ri = sum_{j!=i} |Aij| * s[i] * s[j] を各行で計算
    ri_max = 0.0d0
    
    do i = 1, monoMAT%N
      jS = index(i) + 1
      jE = index(i + 1)
      n1 = monoMAT%n_dof_list(i)
      n2 = monoMAT%n_dof_index(i)

      do dof_i = 1, n1
        ri = 0.0d0
        s_i = s(n2 + dof_i)
        
        do ii = jS, jE
          in = item(ii)
          n3 = monoMAT%n_dof_list(in)
          nz = monoMAT%n_dof_index2(ii)
          aa:do dof_j = 1, n3
            if(n2 + dof_i == monoMAT%n_dof_index(in) + dof_j) cycle aa
            s_j = s(monoMAT%n_dof_index(in) + dof_j)
            Aij_abs = dabs(A(nz + n3*(dof_i-1) + dof_j))
            ri = ri + Aij_abs * s_i * s_j
          enddo aa
        enddo
        
        ri_max = max(ri_max, ri)
      enddo
    enddo

    !# 全プロセス間で ri の最大値を取得
    call monolis_allreduce_R1(ri_max, monolis_mpi_max, monoCOM%comm)

    U = 1.0d0 + ri_max
    omega = 0.9d0 * (2.0d0 / U)

    call monolis_dealloc_R_1d(s)
  end subroutine monolis_solver_JACOBI_get_auto_relax_factor

  subroutine monolis_solver_power_method_eigenvalue_estimation( &
    monoCOM, monoMAT, NNDOF, NPNDOF, lambda_min, lambda_max)
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

    max_iter = 10
    tol = 1.0d-6

    call monolis_alloc_R_1d(v, NPNDOF)
    call monolis_alloc_R_1d(Av, NPNDOF)

    call random_seed()

    do i = 1, NNDOF
      call random_number(v(i))
      v(i) = v(i) - 0.5d0
    enddo

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

    do iter = 1, max_iter
      lambda_old = lambda_max

      ! Av = A * v
      call monolis_matvec_product_main_R(monoCOM, monoMAT, v, Av, tspmv, tcomm)

      lambda = 0.0d0
      norm_v = 0.0d0
      do i = 1, NNDOF
        lambda = lambda + v(i) * Av(i)
        norm_v = norm_v + v(i)**2
      enddo

      call monolis_allreduce_R1(lambda, monolis_mpi_sum, monoCOM%comm)
      call monolis_allreduce_R1(norm_v, monolis_mpi_sum, monoCOM%comm)

      lambda_max = lambda / norm_v

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

      if(iter > 1 .and. dabs(lambda_max - lambda_old) < tol * dabs(lambda_max)) exit
    enddo

    lambda_min = 0.1d0 * lambda_max

    if(lambda_min <= 0.0d0) lambda_min = 0.1d0
    if(lambda_max <= lambda_min) lambda_max = 10.0d0 * lambda_min

    call monolis_dealloc_R_1d(v)
    call monolis_dealloc_R_1d(Av)
  end subroutine monolis_solver_power_method_eigenvalue_estimation

  subroutine monolis_solver_JACOBI_matvec(monoCOM, monoMAT, NNDOF, NPNDOF, X, B, omega, tspmv, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: NNDOF, NPNDOF
    integer(kint) :: i
    real(kdouble) :: X(:), B(:)
    real(kdouble), allocatable :: Y(:)
    real(kdouble), pointer :: D(:)
    real(kdouble) :: omega
    real(kdouble) :: tspmv, tcomm

    D => monoMAT%R%D

    call monolis_alloc_R_1d(Y, NPNDOF)

    call monolis_matvec_product_main_R(monoCOM, monoMAT, X, Y, tspmv, tcomm)

    do i = 1, NNDOF
      X(i) = (1.0d0 - omega)*X(i) + omega*(B(i) - Y(i) + D(i)*X(i)) / D(i)
    enddo

    call monolis_dealloc_R_1d(Y)
  end subroutine monolis_solver_JACOBI_matvec
end module mod_monolis_solver_JACOBI
