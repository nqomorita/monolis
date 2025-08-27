!> ¹ÌÓÐ‚Ž¥½¥ë¥Ð¥â¥¸¥å©`¥ë
module mod_monolis_eigen_solver_wrapper
  use mod_monolis_utils
  use mod_monolis_eigen_solver
  use iso_c_binding

  implicit none

contains

  subroutine monolis_eigen_standard_lanczos_R_c(N, NP, NZ, NDOF, NPNDOF, NZNDOF2, &
    n_dof_list, n_dof_index, n_dof_index2, A, index, item, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    Iarray, Rarray, &
    n_get_eigen, ths, maxiter, val, vec, is_bc) &
    & bind(c, name = "monolis_eigen_standard_lanczos_R_c_main")
    implicit none
    type(monolis_structure) :: monolis
    type(monolis_com) :: monoCOM
    integer(kint_c), intent(in), value :: N, NP, NZ, NDOF
    integer(kint_c), intent(in), value :: NPNDOF
    integer(kint_c), intent(in), value :: NZNDOF2
    integer(kint_c), intent(in), value :: my_rank, comm, comm_size
    integer(kint_c), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(kint_c), intent(in), value :: maxiter
    integer(kint_c), intent(inout), target :: n_get_eigen
    integer(kint_c), intent(in), target :: n_dof_list(NP)
    integer(kint_c), intent(in), target :: n_dof_index(NP + 1)
    integer(kint_c), intent(in), target :: n_dof_index2(NZ + 1)
    integer(kint_c), intent(in), target :: index(NP + 1)
    integer(kint_c), intent(in), target :: item(NZ)
    integer(kint_c), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(kint_c), intent(in), target :: recv_index(recv_n_neib + 1), recv_item(recv_nitem)
    integer(kint_c), intent(in), target :: send_neib_pe(send_n_neib)
    integer(kint_c), intent(in), target :: send_index(send_n_neib + 1), send_item(send_nitem)
    integer(kint_c), intent(in), target :: is_bc(NPNDOF)
    integer(kint_c), intent(in), target :: Iarray(100)
    real(c_double), intent(in), target :: Rarray(100)
    real(c_double), intent(in), target :: A(NZNDOF2)
    real(c_double), intent(in), value :: ths
    real(c_double), intent(out), target :: val(n_get_eigen)
    real(c_double), intent(out), target :: vec(NPNDOF, n_get_eigen)
    integer(kint) :: i
    logical, allocatable :: is_bc_t(:)

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NDOF = NDOF
    monolis%MAT%n_dof_list => n_dof_list
    monolis%MAT%n_dof_index => n_dof_index
    monolis%MAT%n_dof_index2 => n_dof_index2
    monolis%MAT%R%A => A
    monolis%MAT%CSR%index => index
    monolis%MAT%CSR%item => item
    call monolis_palloc_R_1d(monolis%MAT%R%X, NPNDOF)
    call monolis_palloc_R_1d(monolis%MAT%R%B, NPNDOF)

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

    call monolis_alloc_L_1d(is_bc_t, NPNDOF)

    do i = 1, NPNDOF
      is_bc_t(i) = monolis_conv_I2L(is_bc(i))
    enddo

    monoCOM%recv_item = monoCOM%recv_item + 1
    monoCOM%send_item = monoCOM%send_item + 1

    call monolis_eigen_standard_lanczos_R_main( &
      & monolis%PRM, monoCOM, monolis%MAT, n_get_eigen, ths, maxiter, val, vec, is_bc_t)

    monoCOM%recv_item = monoCOM%recv_item - 1
    monoCOM%send_item = monoCOM%send_item - 1
  end subroutine monolis_eigen_standard_lanczos_R_c

  subroutine monolis_eigen_inverted_standard_lanczos_R_c(N, NP, NZ, NDOF, NPNDOF, NZNDOF2, &
    n_dof_list, n_dof_index, n_dof_index2, A, index, item, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    Iarray, Rarray, &
    n_get_eigen, ths, maxiter, val, vec, is_bc) &
    & bind(c, name = "monolis_eigen_inverted_standard_lanczos_R_c_main")
    implicit none
    type(monolis_structure) :: monolis
    type(monolis_com) :: monoCOM
    integer(kint_c), intent(in), value :: N, NP, NZ, NDOF
    integer(kint_c), intent(in), value :: NPNDOF
    integer(kint_c), intent(in), value :: NZNDOF2
    integer(kint_c), intent(in), value :: my_rank, comm, comm_size
    integer(kint_c), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(kint_c), intent(in), value :: maxiter
    integer(kint_c), intent(inout), target :: n_get_eigen
    integer(kint_c), intent(in), target :: n_dof_list(NP)
    integer(kint_c), intent(in), target :: n_dof_index(NP + 1)
    integer(kint_c), intent(in), target :: n_dof_index2(NZ + 1)
    integer(kint_c), intent(in), target :: index(NP + 1)
    integer(kint_c), intent(in), target :: item(NZ)
    integer(kint_c), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(kint_c), intent(in), target :: recv_index(recv_n_neib + 1), recv_item(recv_nitem)
    integer(kint_c), intent(in), target :: send_neib_pe(send_n_neib)
    integer(kint_c), intent(in), target :: send_index(send_n_neib + 1), send_item(send_nitem)
    integer(kint_c), intent(in), target :: is_bc(NPNDOF)
    integer(kint_c), intent(in), target :: Iarray(100)
    real(c_double), intent(in), target :: Rarray(100)
    real(c_double), intent(in), target :: A(NZNDOF2)
    real(c_double), intent(in), value :: ths
    real(c_double), intent(out), target :: val(n_get_eigen)
    real(c_double), intent(out), target :: vec(NPNDOF, n_get_eigen)
    integer(kint) :: i
    logical, allocatable :: is_bc_t(:)

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NDOF = NDOF
    monolis%MAT%n_dof_list => n_dof_list
    monolis%MAT%n_dof_index => n_dof_index
    monolis%MAT%n_dof_index2 => n_dof_index2
    monolis%MAT%R%A => A
    monolis%MAT%CSR%index => index
    monolis%MAT%CSR%item => item
    call monolis_palloc_R_1d(monolis%MAT%R%X, NPNDOF)
    call monolis_palloc_R_1d(monolis%MAT%R%B, NPNDOF)

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

    call monolis_alloc_L_1d(is_bc_t, NPNDOF)

    do i = 1, NPNDOF
      is_bc_t(i) = monolis_conv_I2L(is_bc(i))
    enddo

    monoCOM%recv_item = monoCOM%recv_item + 1
    monoCOM%send_item = monoCOM%send_item + 1

    call monolis_eigen_inverted_standard_lanczos_R_main( &
      & monolis%PRM, monoCOM, monolis%MAT, monolis%PREC, n_get_eigen, ths, maxiter, val, vec, is_bc_t)

    monoCOM%recv_item = monoCOM%recv_item - 1
    monoCOM%send_item = monoCOM%send_item - 1
  end subroutine monolis_eigen_inverted_standard_lanczos_R_c

  subroutine monolis_get_condition_number_R_c(N, NP, NZ, NDOF, NPNDOF, NZNDOF2, &
    n_dof_list, n_dof_index, n_dof_index2, A, index, item, &
    my_rank, comm, comm_size, n_internal_vertex, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    Iarray, Rarray, &
    singular_value_max, singular_value_min) &
    & bind(c, name = "monolis_get_condition_number_R_c_main")
    implicit none
    type(monolis_structure) :: monolis
    type(monolis_com) :: monoCOM
    integer(kint_c), intent(in), value :: N, NP, NZ, NDOF
    integer(kint_c), intent(in), value :: NPNDOF
    integer(kint_c), intent(in), value :: NZNDOF2
    integer(kint_c), intent(in), value :: my_rank, comm, comm_size, n_internal_vertex
    integer(kint_c), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(kint_c), intent(in), target :: n_dof_list(NP)
    integer(kint_c), intent(in), target :: n_dof_index(NP + 1)
    integer(kint_c), intent(in), target :: n_dof_index2(NZ + 1)
    integer(kint_c), intent(in), target :: index(NP + 1)
    integer(kint_c), intent(in), target :: item(NZ)
    integer(kint_c), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(kint_c), intent(in), target :: recv_index(recv_n_neib + 1), recv_item(recv_nitem)
    integer(kint_c), intent(in), target :: send_neib_pe(send_n_neib)
    integer(kint_c), intent(in), target :: send_index(send_n_neib + 1), send_item(send_nitem)
    integer(kint_c), intent(in), target :: Iarray(100)
    real(c_double), intent(in), target :: Rarray(100)
    real(c_double), intent(in), target :: A(NZNDOF2)
    real(c_double), intent(out), target :: singular_value_max
    real(c_double), intent(out), target :: singular_value_min
    integer(kint) :: i

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NDOF = NDOF
    monolis%MAT%n_dof_list => n_dof_list
    monolis%MAT%n_dof_index => n_dof_index
    monolis%MAT%n_dof_index2 => n_dof_index2
    monolis%MAT%R%A => A
    monolis%MAT%CSR%index => index
    monolis%MAT%CSR%item => item
    call monolis_palloc_R_1d(monolis%MAT%R%X, NPNDOF)
    call monolis_palloc_R_1d(monolis%MAT%R%B, NPNDOF)

    !> for monoCOM
    monoCOM%my_rank = my_rank
    monoCOM%comm = comm
    monoCOM%comm_size = comm_size
    monoCOM%n_internal_vertex = n_internal_vertex
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

    call monolis_get_condition_number_R(monolis, monoCOM, singular_value_max, singular_value_min)

    monoCOM%recv_item = monoCOM%recv_item - 1
    monoCOM%send_item = monoCOM%send_item - 1
  end subroutine monolis_get_condition_number_R_c
end module mod_monolis_eigen_solver_wrapper
