!> ���Υ���Х⥸��`��
module mod_monolis_solve_wrapper
  use mod_monolis_utils
  use mod_monolis_solve
  use iso_c_binding

  implicit none

contains

  subroutine monolis_solve_R_c(N, NP, NZ, NDOF, A, X, B, index, item, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    Iarray, Rarray) &
    & bind(c, name = "monolis_solve_R_c_main")
    implicit none
    type(monolis_structure) :: monolis
    integer(c_int), intent(in), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), value :: my_rank, comm, comm_size
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), target :: index(NP + 1)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(recv_n_neib + 1), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(send_n_neib + 1), send_item(send_nitem)
    integer(c_int), intent(in), target :: Iarray(100)
    real(c_double), intent(in), target :: Rarray(100)
    real(c_double), intent(in), target :: A(NDOF*NDOF*NZ)
    real(c_double), intent(in), target :: X(NDOF*NP)
    real(c_double), intent(in), target :: B(NDOF*NP)

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NDOF = NDOF
    !monolis%MAT%A => A
    !monolis%MAT%X => X
    !monolis%MAT%B => B
    !monolis%MAT%index => index
    !monolis%MAT%item => item

    !> for monoCOM
    monoliS%COM%my_rank = my_rank
    monoliS%COM%comm = comm
    monoliS%COM%comm_size = comm_size
    monoliS%COM%recv_n_neib = recv_n_neib
    !monoliS%COM%recv_neib_pe => recv_neib_pe
    !monoliS%COM%recv_index => recv_index
    !monoliS%COM%recv_item => recv_item
    !monoliS%COM%send_n_neib = send_n_neib
    !monoliS%COM%send_neib_pe => send_neib_pe
    !monoliS%COM%send_index => send_index
    !monoliS%COM%send_item => send_item

    call monolis_solve_(monolis%PRM, monolis%COM, monolis%MAT, monolis%PREC)
  end subroutine monolis_solve_R_c

end module mod_monolis_solve_wrapper