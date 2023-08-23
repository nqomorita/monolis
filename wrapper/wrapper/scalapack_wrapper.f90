!> ���Υ���Х⥸��`��
module mod_monolis_scalapack_wrapper
  use mod_monolis_utils
  use mod_monolis_scalapack
  use iso_c_binding

  implicit none

contains

  subroutine monolis_scalapack_gesvd_R_c(N_loc, M, P, A, S, V, D, comm) &
    & bind(c, name = "monolis_scalapack_gesvd_R_c_main")
    implicit none
    !> ���Фδ󤭤������� N��
    integer(c_int), value :: N_loc
    !> ���Фδ󤭤������� M��
    integer(c_int), value :: M
    !> ���Фδ󤭤�����С��
    integer(c_int), value :: P
    !> �������У�N_loc x M��
    real(c_double), target :: A(N_loc*M)
    !> ���خ����У�N_loc x P��
    real(c_double), target :: S(N_loc*P)
    !> �خ�����P��
    real(c_double), target :: V(P)
    !> ���خ����У�P x M��
    real(c_double), target :: D(P*M)
    !> ���ߥ�˥��`��
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
