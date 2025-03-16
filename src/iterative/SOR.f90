!> SOR 法モジュール
module mod_monolis_solver_SOR
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
  public monolis_solver_SOR

contains

  subroutine monolis_solver_SOR(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), target, intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: N, NP, NDOF
    integer(kint) :: iter
    real(kdouble) :: R2, B2
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    real(kdouble), pointer :: B(:), X(:)
    real(kdouble), allocatable :: R(:)
    logical :: is_converge

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF

    X => monoMAT%R%X
    B => monoMAT%R%B

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_alloc_R_1d(R, NDOF*NP)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    call monolis_solver_SOR_setup(monoMAT, monoPREC)
    call monolis_inner_product_main_R(monoCOM, N, NDOF, B, B, B2, tdotp, tcomm_dotp)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_solver_SOR_matvec(monoCOM, monoMAT, NDOF, X, B, tspmv, tcomm_spmv)
      call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R, R, R2, tdotp, tcomm_dotp)
      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit
    enddo

    call monolis_mpi_update_R(monoCOM, NDOF, X, tcomm_spmv)

    call monolis_dealloc_R_1d(R)
  end subroutine monolis_solver_SOR

  subroutine monolis_solver_SOR_setup(monoMAT, monoPREC)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoPREC
    integer(kint) :: i, ii, j, jS, jE, in, k, l, N, NP, NDOF, NDOF2
    integer(kint), pointer :: index(:), item(:)
    real(kdouble), pointer :: A(:), ALU(:)
    real(kdouble), allocatable :: T(:), LU(:,:)

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    A => monoMAT%R%A
    index => monoMAT%CSR%index
    item  => monoMAT%CSR%item

    !allocate(T(NDOF))
    !allocate(LU(NDOF,NDOF))
    !allocate(monoMAT%monoTree%D(NDOF2*NP))
    !ALU => monoMAT%monoTree%D
    !T   = 0.0d0
    !ALU = 0.0d0
    !LU  = 0.0d0

    do i = 1, N
      jS = index(i-1) + 1
      jE = index(i)
      do ii = jS, jE
        in = item(ii)
        if(i == in)then
          do j = 1, NDOF
            do k = 1, NDOF
              LU(j,k) = A(NDOF2*(i-1) + NDOF*(j-1) + k)
            enddo
          enddo
          do k = 1, NDOF
            LU(k,k) = 1.0d0/LU(k,k)
            do l = k+1, NDOF
              LU(l,k) = LU(l,k)*LU(k,k)
              do j = k+1, NDOF
                T(j) = LU(l,j) - LU(l,k)*LU(k,j)
              enddo
              do j = k+1, NDOF
                LU(l,j) = T(j)
              enddo
            enddo
          enddo
          do j = 1, NDOF
            do k = 1, NDOF
              ALU(NDOF2*(i-1) + NDOF*(j-1) + k) = LU(j,k)
            enddo
          enddo
        endif
      enddo
    enddo

    !deallocate(T)
    !deallocate(LU)
  end subroutine monolis_solver_SOR_setup

  subroutine monolis_solver_SOR_matvec(monoCOM, monoMAT, NDOF, X, B, tspmv, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, k, l, in, N, NDOF, NDOF2, jS, jE
    integer(kint), pointer :: index(:), item(:)
    real(kdouble) :: X(:), B(:), XT(NDOF), YT(NDOF), DT(NDOF), WT(NDOF)
    real(kdouble), pointer :: A(:), ALU(:)
    real(kdouble) :: t1, t2, omega
    real(kdouble) :: tspmv, tcomm

    !ALU => monoMAT%monoTree%D
    !N     = monoMAT%N
    !NDOF  = monoMAT%NDOF
    !NDOF2 = NDOF*NDOF
    !A => monoMAT%R%A
    !index => monoMAT%index
    !item  => monoMAT%item
    !omega = 1.0d0

    call monolis_mpi_update_R(monoCOM, NDOF, X, tcomm)

    do i = 1, N
      DT = 0.0d0
      do k = 1, NDOF
        XT(k) = X(NDOF*(i-1) + k)
      enddo
      do j = 1, NDOF
        do k = 1, NDOF
          DT(j) = DT(j) + A(NDOF2*(i-1) + NDOF*(j-1) + k)*XT(k)
        enddo
      enddo

      YT = 0.0d0
      jS = index(i-1) + 1
      jE = index(i  )
      do j = jS, jE
        in = item(j)
        if(in < i)then
          do k = 1, NDOF
            XT(k) = X(NDOF*(in-1) + k)
          enddo
          do k = 1, NDOF
            do l = 1, NDOF
              YT(k) = YT(k) - A(NDOF2*(j-1) + NDOF*(k-1) + l)*XT(l)
            enddo
          enddo
        endif
      enddo

      jS = index(i-1) + 1
      jE = index(i  )
      do j = jS, jE
        in = item(j)
        if(i < in)then
          do k = 1, NDOF
            XT(k) = X(NDOF*(in-1) + k)
          enddo
          do k = 1, NDOF
            do l = 1, NDOF
              YT(k) = YT(k) - A(NDOF2*(j-1) + NDOF*(k-1) + l)*XT(l)
            enddo
          enddo
        endif
      enddo

      do k = 1, NDOF
        WT(k) = omega*B(NDOF*(i-1) + k) + (1.0d0 - omega)*DT(k) + omega*YT(k)
      enddo

      do j = 2, NDOF
        do k = 1, j-1
          WT(j) = WT(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*WT(k)
        enddo
      enddo
      do j = NDOF, 1, -1
        do k = NDOF, j+1, -1
          WT(j) = WT(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*WT(k)
        enddo
        WT(j) = ALU(NDOF2*(i-1) + (NDOF+1)*(j-1) + 1)*WT(j)
      enddo
      do k = 1, NDOF
        X(NDOF*(i-1) + k) = WT(k)
      enddo
    enddo
  end subroutine monolis_solver_SOR_matvec
end module mod_monolis_solver_SOR
