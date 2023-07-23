!> ¹ÌÓÐ‚Ž¥½¥ë¥Ð¥â¥¸¥å©`¥ë
module mod_monolis_eigen_solver_wrapper
  use mod_monolis_utils
  use mod_monolis_eigen_solver
  use iso_c_binding

  implicit none

contains

  subroutine monolis_eigen_standard_lanczos_R_c(N, NP, NZ, NDOF, A, index, item, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    Iarray, Rarray, &
    n_get_eigen, ths, maxiter, val, vec, is_bc) &
    & bind(c, name = "monolis_eigen_standard_lanczos_R_c_main")
    implicit none
    type(monolis_structure) :: monolis
    type(monolis_com) :: monoCOM
    integer(c_int), intent(in), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), value :: my_rank, comm, comm_size
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), value :: maxiter
    integer(c_int), intent(inout), target :: n_get_eigen
    integer(c_int), intent(in), target :: index(NP + 1)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(recv_n_neib + 1), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(send_n_neib + 1), send_item(send_nitem)
    integer(c_int), intent(in), target :: is_bc(NDOF*NP)
    integer(c_int), intent(in), target :: Iarray(100)
    real(c_double), intent(in), target :: Rarray(100)
    real(c_double), intent(in), target :: A(NDOF*NDOF*NZ)
    real(c_double), intent(in), value :: ths
    real(c_double), intent(out), target :: val(n_get_eigen)
    real(c_double), intent(out), target :: vec(NDOF*NP, n_get_eigen)
    integer(kint) :: i
    logical, allocatable :: is_bc_t(:)

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NDOF = NDOF
    monolis%MAT%R%A => A
    monolis%MAT%CSR%index => index
    monolis%MAT%CSR%item => item
    call monolis_palloc_R_1d(monolis%MAT%R%X, NP*NDOF)
    call monolis_palloc_R_1d(monolis%MAT%R%B, NP*NDOF)

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

    call monolis_alloc_L_1d(is_bc_t, NP*NDOF)

    do i = 1, NP*NDOF
      is_bc_t(i) = monolis_conv_I2L(is_bc(i))
    enddo

    monoCOM%recv_item = monoCOM%recv_item + 1
    monoCOM%send_item = monoCOM%send_item + 1

    call monolis_eigen_standard_lanczos_R_main( &
      & monolis%PRM, monoCOM, monolis%MAT, n_get_eigen, ths, maxiter, val, vec, is_bc_t)

    monoCOM%recv_item = monoCOM%recv_item - 1
    monoCOM%send_item = monoCOM%send_item - 1
  end subroutine monolis_eigen_standard_lanczos_R_c

  subroutine monolis_eigen_inverted_standard_lanczos_R_c(N, NP, NZ, NDOF, A, index, item, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    Iarray, Rarray, &
    n_get_eigen, ths, maxiter, val, vec, is_bc) &
    & bind(c, name = "monolis_eigen_inverted_standard_lanczos_R_c_main")
    implicit none
    type(monolis_structure) :: monolis
    type(monolis_com) :: monoCOM
    integer(c_int), intent(in), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), value :: my_rank, comm, comm_size
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), value :: maxiter
    integer(c_int), intent(inout), target :: n_get_eigen
    integer(c_int), intent(in), target :: index(NP + 1)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(recv_n_neib + 1), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(send_n_neib + 1), send_item(send_nitem)
    integer(c_int), intent(in), target :: is_bc(NDOF*NP)
    integer(c_int), intent(in), target :: Iarray(100)
    real(c_double), intent(in), target :: Rarray(100)
    real(c_double), intent(in), target :: A(NDOF*NDOF*NZ)
    real(c_double), intent(in), value :: ths
    real(c_double), intent(out), target :: val(n_get_eigen)
    real(c_double), intent(out), target :: vec(NDOF*NP, n_get_eigen)
    integer(kint) :: i
    logical, allocatable :: is_bc_t(:)

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NDOF = NDOF
    monolis%MAT%R%A => A
    monolis%MAT%CSR%index => index
    monolis%MAT%CSR%item => item
    call monolis_palloc_R_1d(monolis%MAT%R%X, NP*NDOF)
    call monolis_palloc_R_1d(monolis%MAT%R%B, NP*NDOF)

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

    call monolis_alloc_L_1d(is_bc_t, NP*NDOF)

    do i = 1, NP*NDOF
      is_bc_t(i) = monolis_conv_I2L(is_bc(i))
    enddo

    monoCOM%recv_item = monoCOM%recv_item + 1
    monoCOM%send_item = monoCOM%send_item + 1

    call monolis_eigen_inverted_standard_lanczos_R_main( &
      & monolis%PRM, monoCOM, monolis%MAT, monolis%PREC, n_get_eigen, ths, maxiter, val, vec, is_bc_t)

    monoCOM%recv_item = monoCOM%recv_item - 1
    monoCOM%send_item = monoCOM%send_item - 1
  end subroutine monolis_eigen_inverted_standard_lanczos_R_c
end module mod_monolis_eigen_solver_wrapper
