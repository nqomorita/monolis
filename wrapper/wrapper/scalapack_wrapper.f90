!> 線形ソルバモジュール
module mod_monolis_scalapack_wrapper
  use mod_monolis_utils
  use mod_monolis_scalapack
  use iso_c_binding

  implicit none

contains

  subroutine monolis_scalapack_gesvd_R_c(N_loc, M, P, A, S, V, D, comm) &
    & bind(c, name = "monolis_scalapack_gesvd_R_c_main")
    implicit none
    !> 行列の大きさ（行数 N）
    integer(c_int), value :: N_loc
    !> 行列の大きさ（列数 M）
    integer(c_int), value :: M
    !> 行列の大きさの最小値
    integer(c_int), value :: P
    !> 入力行列（N_loc x M）
    real(c_double), target :: A(N_loc*M)
    !> 左特異行列（N_loc x P）
    real(c_double), target :: S(N_loc*P)
    !> 特異値（P）
    real(c_double), target :: V(P)
    !> 右特異行列（P x M）
    real(c_double), target :: D(P*M)
    !> コミュニケータ
    integer(c_int), value :: comm
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

    call monolis_scalapack_gesvd_R(N_loc, M, A_temp, S_temp, V, D_temp, comm)

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
end module mod_monolis_scalapack_wrapper
