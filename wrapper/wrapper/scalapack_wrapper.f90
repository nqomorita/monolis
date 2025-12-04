!> 線形ソルバモジュール
module mod_monolis_scalapack_wrapper
  use mod_monolis_utils
  use mod_monolis_scalapack
  use iso_c_binding

  implicit none

contains

  subroutine monolis_scalapack_comm_initialize_c(comm, scalapack_comm) &
    & bind(c, name = "monolis_scalapack_comm_initialize_c_main")
    implicit none
    !> [in] コミュニケータ
    integer(kint_c), value :: comm
    !> [in] コミュニケータ
    integer(kint_c) :: scalapack_comm
    call monolis_scalapack_comm_initialize(comm, scalapack_comm)
  end subroutine monolis_scalapack_comm_initialize_c

  subroutine monolis_scalapack_comm_finalize_c(scalapack_comm) &
    & bind(c, name = "monolis_scalapack_comm_finalize_c_main")
    implicit none
    !> [in] コミュニケータ
    integer(kint_c), value :: scalapack_comm
    call monolis_scalapack_comm_finalize(scalapack_comm)
  end subroutine monolis_scalapack_comm_finalize_c

  subroutine monolis_scalapack_gesvd_R_c(N_loc, M, P, A, S, V, D, comm, scalapack_comm) &
    & bind(c, name = "monolis_scalapack_gesvd_R_c_main")
    implicit none
    !> 行列の大きさ（行数 N）
    integer(kint_c), value :: N_loc
    !> 行列の大きさ（列数 M）
    integer(kint_c), value :: M
    !> 行列の大きさの最小値
    integer(kint_c), value :: P
    !> 入力行列（N_loc x M）
    real(c_double), target :: A(N_loc*M)
    !> 左特異行列（N_loc x P）
    real(c_double), target :: S(N_loc*P)
    !> 特異値（P）
    real(c_double), target :: V(P)
    !> 右特異行列（P x M）
    real(c_double), target :: D(P*M)
    !> コミュニケータ
    integer(kint_c), value :: comm
    !> scalapack コミュニケータ
    integer(kint_c), value :: scalapack_comm
    integer(kint) :: i, j
    real(kdouble), allocatable :: A_temp(:,:)
    real(kdouble), allocatable :: S_temp(:,:)
    real(kdouble), allocatable :: D_temp(:,:)

    call monolis_alloc_R_2d(A_temp, N_loc, M)
    call monolis_alloc_R_2d(S_temp, N_loc, P)
    call monolis_alloc_R_2d(D_temp, P, M)

    do i = 1, M
      do j = 1, N_loc
        A_temp(j,i) = A(N_loc*(i-1) + j)
      enddo
    enddo

    call monolis_scalapack_gesvd_R(N_loc, M, A_temp, S_temp, V, D_temp, comm, scalapack_comm)

    do i = 1, P
      do j = 1, N_loc
        S(N_loc*(i-1) + j) = S_temp(j,i)
      enddo
    enddo

    do i = 1, M
      do j = 1, P
         D(P*(i-1) + j) = D_temp(j,i)
      enddo
    enddo
  end subroutine monolis_scalapack_gesvd_R_c

  subroutine monolis_scalapack_getrf_R_c(N_loc, N, A, ipiv, comm, scalapack_comm) &
    & bind(c, name = "monolis_scalapack_getrf_R_c_main")
    implicit none
    !> 行列の大きさ（ローカル行数）
    integer(kint_c), value :: N_loc
    !> 行列の大きさ（全体のサイズ N x N）
    integer(kint_c), value :: N
    !> 入力行列（N_loc x N）
    real(c_double), target :: A(N_loc*N)
    !> ピボット情報（N_loc）
    integer(kint_c), target :: ipiv(N_loc)
    !> コミュニケータ
    integer(kint_c), value :: comm
    !> scalapack コミュニケータ
    integer(kint_c), value :: scalapack_comm
    integer(kint) :: i, j
    real(kdouble), allocatable :: A_temp(:,:)
    integer(kint), allocatable :: ipiv_temp(:)

    call monolis_alloc_R_2d(A_temp, N_loc, N)
    call monolis_alloc_I_1d(ipiv_temp, N_loc)

    do i = 1, N
      do j = 1, N_loc
        A_temp(j,i) = A(N_loc*(i-1) + j)
      enddo
    enddo

    call monolis_scalapack_getrf_R(N_loc, N, A_temp, ipiv_temp, comm, scalapack_comm)

    do i = 1, N
      do j = 1, N_loc
        A(N_loc*(i-1) + j) = A_temp(j,i)
      enddo
    enddo

    do i = 1, N_loc
      ipiv(i) = ipiv_temp(i)
    enddo
  end subroutine monolis_scalapack_getrf_R_c

  subroutine monolis_scalapack_getrs_R_c(N_loc, N, NRHS, A, ipiv, B, comm, scalapack_comm) &
    & bind(c, name = "monolis_scalapack_getrs_R_c_main")
    implicit none
    !> 行列の大きさ（ローカル行数）
    integer(kint_c), value :: N_loc
    !> 行列の大きさ（全体のサイズ N x N）
    integer(kint_c), value :: N
    !> 右辺ベクトルの数
    integer(kint_c), value :: NRHS
    !> LU分解された行列（N_loc x N）
    real(c_double), target :: A(N_loc*N)
    !> ピボット情報（N_loc）
    integer(kint_c), target :: ipiv(N_loc)
    !> 右辺ベクトル（N_loc x NRHS）
    real(c_double), target :: B(N_loc*NRHS)
    !> コミュニケータ
    integer(kint_c), value :: comm
    !> scalapack コミュニケータ
    integer(kint_c), value :: scalapack_comm
    integer(kint) :: i, j
    real(kdouble), allocatable :: A_temp(:,:)
    real(kdouble), allocatable :: B_temp(:,:)
    integer(kint), allocatable :: ipiv_temp(:)

    call monolis_alloc_R_2d(A_temp, N_loc, N)
    call monolis_alloc_R_2d(B_temp, N_loc, NRHS)
    call monolis_alloc_I_1d(ipiv_temp, N_loc)

    do i = 1, N
      do j = 1, N_loc
        A_temp(j,i) = A(N_loc*(i-1) + j)
      enddo
    enddo

    do i = 1, NRHS
      do j = 1, N_loc
        B_temp(j,i) = B(N_loc*(i-1) + j)
      enddo
    enddo

    do i = 1, N_loc
      ipiv_temp(i) = ipiv(i)
    enddo

    call monolis_scalapack_getrs_R(N_loc, N, NRHS, A_temp, ipiv_temp, B_temp, comm, scalapack_comm)

    do i = 1, NRHS
      do j = 1, N_loc
        B(N_loc*(i-1) + j) = B_temp(j,i)
      enddo
    enddo
  end subroutine monolis_scalapack_getrs_R_c
end module mod_monolis_scalapack_wrapper
