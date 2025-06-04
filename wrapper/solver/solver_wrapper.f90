!> ���Υ���Х⥸��`��
module mod_monolis_solve_wrapper
  use mod_monolis_utils
  use mod_monolis_solve
  use iso_c_binding

  implicit none

contains

  subroutine monolis_solve_R_c(N, NP, NZ, NDOF, NPNDOF, NZNDOF2, &
    n_dof_list, n_dof_index, n_dof_index2, &
    A, X, B, index, item, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    Iarray, Rarray) &
    & bind(c, name = "monolis_solve_R_c_main")
    implicit none
    type(monolis_structure) :: monolis
    type(monolis_com) :: monoCOM
    integer(c_int), intent(in), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), value :: NPNDOF
    integer(c_int), intent(in), value :: NZNDOF2
    integer(c_int), intent(in), value :: my_rank, comm, comm_size
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), target :: n_dof_list(NP)
    integer(c_int), intent(in), target :: n_dof_index(NP + 1)
    integer(c_int), intent(in), target :: n_dof_index2(NZ + 1)
    integer(c_int), intent(in), target :: index(NP + 1)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(recv_n_neib + 1), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(send_n_neib + 1), send_item(send_nitem)
    integer(c_int), target :: Iarray(100)
    real(c_double), target :: Rarray(100)
    real(c_double), intent(in), target :: A(NZNDOF2)
    real(c_double), target :: X(NPNDOF)
    real(c_double), intent(in), target :: B(NPNDOF)
    integer(kint) :: i

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NDOF = NDOF
    monolis%MAT%n_dof_list => n_dof_list
    monolis%MAT%n_dof_index => n_dof_index
    monolis%MAT%n_dof_index2 => n_dof_index2
    monolis%MAT%R%A => A
    monolis%MAT%R%X => X
    monolis%MAT%R%B => B
    monolis%MAT%CSR%index => index
    monolis%MAT%CSR%item => item

    !> for monoCOM
    monoCOM%my_rank = my_rank
    monoCOM%comm = comm
    monoCOM%comm_size = comm_size
    monoCOM%recv_n_neib = recv_n_neib
    monoCOM%recv_neib_pe => recv_neib_pe
    monoCOM%recv_index => recv_index
    monoCOM%recv_item => recv_item
    monoCOM%send_n_neib = send_n_neib
    monoCOM%send_neib_pe => send_neib_pe
    monoCOM%send_index => send_index
    monoCOM%send_item => send_item

    do i = 1, monolis_prm_Iarray_size - 1
      monoliS%PRM%Iarray(i) = Iarray(i + 1)
    enddo

    do i = 1, monolis_prm_Rarray_size - 1
      monoliS%PRM%Rarray(i) = Rarray(i + 1)
    enddo

    monoCOM%recv_item = monoCOM%recv_item + 1
    monoCOM%send_item = monoCOM%send_item + 1

    call monolis_solve_main_R(monolis%PRM, monoCOM, monolis%MAT, monolis%PREC)

    monoCOM%recv_item = monoCOM%recv_item - 1
    monoCOM%send_item = monoCOM%send_item - 1

    do i = 1, monolis_prm_Iarray_size - 1
      Iarray(i + 1) = monoliS%PRM%Iarray(i)
    enddo

    do i = 1, monolis_prm_Rarray_size - 1
      Rarray(i + 1) = monoliS%PRM%Rarray(i)
    enddo
  end subroutine monolis_solve_R_c

  subroutine monolis_solve_C_c(N, NP, NZ, NDOF, NPNDOF, NZNDOF2, &
    n_dof_list, n_dof_index, n_dof_index2, &
    A, X, B, index, item, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    Iarray, Rarray) &
    & bind(c, name = "monolis_solve_C_c_main")
    implicit none
    type(monolis_structure) :: monolis
    type(monolis_com) :: monoCOM
    integer(c_int), intent(in), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), value :: NPNDOF
    integer(c_int), intent(in), value :: NZNDOF2
    integer(c_int), intent(in), value :: my_rank, comm, comm_size
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), target :: n_dof_list(NP)
    integer(c_int), intent(in), target :: n_dof_index(NP + 1)
    integer(c_int), intent(in), target :: n_dof_index2(NZ + 1)
    integer(c_int), intent(in), target :: index(NP + 1)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(recv_n_neib + 1), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(send_n_neib + 1), send_item(send_nitem)
    integer(c_int), target :: Iarray(100)
    real(c_double), target :: Rarray(100)
    complex(c_double), intent(in), target :: A(NZNDOF2)
    complex(c_double), target :: X(NPNDOF)
    complex(c_double), intent(in), target :: B(NPNDOF)
    integer(kint) :: i

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NDOF = NDOF
    monolis%MAT%n_dof_list => n_dof_list
    monolis%MAT%n_dof_index => n_dof_index
    monolis%MAT%n_dof_index2 => n_dof_index2
    monolis%MAT%C%A => A
    monolis%MAT%C%X => X
    monolis%MAT%C%B => B
    monolis%MAT%CSR%index => index
    monolis%MAT%CSR%item => item

    !> for monoCOM
    monoCOM%my_rank = my_rank
    monoCOM%comm = comm
    monoCOM%comm_size = comm_size
    monoCOM%recv_n_neib = recv_n_neib
    monoCOM%recv_neib_pe => recv_neib_pe
    monoCOM%recv_index => recv_index
    monoCOM%recv_item => recv_item
    monoCOM%send_n_neib = send_n_neib
    monoCOM%send_neib_pe => send_neib_pe
    monoCOM%send_index => send_index
    monoCOM%send_item => send_item

    do i = 1, monolis_prm_Iarray_size - 1
      monoliS%PRM%Iarray(i) = Iarray(i + 1)
    enddo

    do i = 1, monolis_prm_Rarray_size - 1
      monoliS%PRM%Rarray(i) = Rarray(i + 1)
    enddo

    monoCOM%recv_item = monoCOM%recv_item + 1
    monoCOM%send_item = monoCOM%send_item + 1

    call monolis_solve_main_C(monolis%PRM, monoCOM, monolis%MAT, monolis%PREC)

    monoCOM%recv_item = monoCOM%recv_item - 1
    monoCOM%send_item = monoCOM%send_item - 1

    do i = 1, monolis_prm_Iarray_size - 1
      Iarray(i + 1) = monoliS%PRM%Iarray(i)
    enddo

    do i = 1, monolis_prm_Rarray_size - 1
      Rarray(i + 1) = monoliS%PRM%Rarray(i)
    enddo
  end subroutine monolis_solve_C_c
end module mod_monolis_solve_wrapper
