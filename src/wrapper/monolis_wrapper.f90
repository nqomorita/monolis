module mod_monolis_wrapper
  use iso_c_binding
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_sparse_util
  use mod_monolis_stdlib
  implicit none

contains

  !> initialize section
  subroutine monolis_global_initialize_c() &
    & bind(c, name = "monolis_global_initialize")
    implicit none
    call monolis_global_initialize()
  end subroutine monolis_global_initialize_c

  subroutine monolis_global_finalize_c() &
    & bind(c, name = "monolis_global_finalize")
    implicit none
    call monolis_global_finalize()
  end subroutine monolis_global_finalize_c

  !> getter
  function monolis_get_global_comm_c()&
    & bind(c, name = "monolis_get_global_comm")
    implicit none
    integer(c_int) :: monolis_get_global_comm_c
    monolis_get_global_comm_c =  monolis_global_comm()
  end function monolis_get_global_comm_c

  function monolis_get_global_commsize_c()&
    & bind(c, name = "monolis_get_global_commsize")
    implicit none
    integer(c_int) :: monolis_get_global_commsize_c
    monolis_get_global_commsize_c =  monolis_global_commsize()
  end function monolis_get_global_commsize_c

  function monolis_get_global_myrank_c()&
    & bind(c, name = "monolis_get_global_myrank")
    implicit none
    integer(c_int) :: monolis_get_global_myrank_c
    monolis_get_global_myrank_c =  monolis_global_myrank()
  end function monolis_get_global_myrank_c

  function monolis_get_time_c()&
    & bind(c, name = "monolis_get_time")
    implicit none
    real(c_double) :: monolis_get_time_c
    monolis_get_time_c =  monolis_get_time()
  end function monolis_get_time_c

  function monolis_get_time_sync_c()&
    & bind(c, name = "monolis_get_time_sync")
    implicit none
    real(c_double) :: monolis_get_time_sync_c
    monolis_get_time_sync_c =  monolis_get_time_sync()
  end function monolis_get_time_sync_c

  subroutine monolis_barrier_c(comm)&
    & bind(c, name = "monolis_barrier_c_main")
    use mod_monolis_com
    implicit none
    integer(c_int), intent(in), value :: comm
    call monolis_barrier_(comm)
  end subroutine monolis_barrier_c

  !> mat
  subroutine monolis_get_CRR_format_c(N, NZ, index, item, indexR, itemR, permR) &
    & bind(c, name = "monolis_get_CRR_format")
    implicit none
    integer(c_int), intent(in), value :: N, NZ
    integer(c_int), intent(in), target :: index(0:N)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), target :: indexR(0:N)
    integer(c_int), target :: itemR(NZ)
    integer(c_int), target :: permR(NZ)
    integer(kint), pointer :: indexRt(:)
    integer(kint), pointer :: itemRt(:)
    integer(kint), pointer :: permRt(:)

    indexRt => indexR
    itemRt => itemR
    permRt => permR
    call monolis_get_CRR_format(N, NZ, index, item, indexRt, itemRt, permRt)
  end subroutine monolis_get_CRR_format_c

  subroutine monolis_sparse_matrix_add_value_c(N, NZ, NDOF, index, item, A, i, j, sub_i, sub_j, val) &
    & bind(c, name = "monolis_add_scalar_to_sparse_matrix_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, i, j, sub_i, sub_j
    integer(c_int), intent(in), target :: index(0:N)
    integer(c_int), intent(in), target :: item(NZ)
    real(c_double), target :: A(NDOF*NDOF*NZ)
    real(c_double), value :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_sparse_matrix_add_value(index, item, A, NDOF, i_t, j_t, sub_i_t, sub_j_t, val)
  end subroutine monolis_sparse_matrix_add_value_c

  subroutine monolis_add_sparse_matrix_c(N, NZ, NDOF, NBF, index, item, A, conn, mat) &
    & bind(c, name = "monolis_add_sparse_matrix_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, NBF
    integer(c_int), intent(in), target :: index(0:N)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: conn(NBF)
    real(c_double), target :: A(NDOF*NDOF*NZ)
    real(c_double), target :: mat(NDOF*NDOF*NBF*NBF)
    integer(kint) :: conn_t(NBF)

    conn_t = conn + 1
    call monolis_sparse_matrix_add_matrix(index, item, A, NBF, NBF, NDOF, conn_t, conn_t, mat)
  end subroutine monolis_add_sparse_matrix_c

  subroutine monolis_get_scalar_from_sparse_matrix_c(N, NZ, NDOF, index, item, A, i, j, sub_i, sub_j, val, is_find) &
    & bind(c, name = "monolis_get_scalar_from_sparse_matrix_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, i, j, sub_i, sub_j
    integer(c_int), intent(in), target :: index(0:N)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), target :: is_find
    real(c_double), target :: A(NDOF*NDOF*NZ)
    real(c_double), target :: val
    integer(kint) :: i_t, j_t, sub_i_t, sub_j_t
    logical :: is_find_t

    i_t = i + 1
    j_t = j + 1
    sub_i_t = sub_i + 1
    sub_j_t = sub_j + 1
    call monolis_sparse_matrix_get_value(index, item, A, NDOF, i_t, j_t, sub_i_t, sub_j_t, val, is_find_t)
    is_find = 0
    if(is_find_t) is_find = 1
  end subroutine monolis_get_scalar_from_sparse_matrix_c

  subroutine monolis_set_Dirichlet_bc_c(N, NZ, NDOF, index, item, indexR, itemR, permR, &
    & A, B, nid, ndof_bc, val) &
    & bind(c, name = "monolis_set_Dirichlet_bc_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, nid, ndof_bc
    integer(c_int), intent(in), target :: index(0:N)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: indexR(0:N)
    integer(c_int), intent(in), target :: itemR(NZ)
    integer(c_int), intent(in), target :: permR(NZ)
    real(c_double), intent(in), value :: val
    real(c_double), target :: A(NDOF*NDOF*NZ)
    real(c_double), target :: B(NDOF*N)
    integer(kint) :: nid_t, ndof_bc_t

    nid_t = nid + 1
    ndof_bc_t = ndof_bc + 1
    call monolis_sparse_matrix_add_bc(index, item, A, B, indexR, itemR, permR, &
      & ndof, nid_t, ndof_bc_t, val)
  end subroutine monolis_set_Dirichlet_bc_c

  !> linenar alg
  subroutine monolis_matvec_product_c(N, NP, NZ, NDOF, A, X, Y, index, item, &
    myrank, comm, commsize, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item) &
    & bind(c, name = "monolis_matvec_product_c_main")
    use mod_monolis_com
    use mod_monolis_matvec
    implicit none
    type(monolis_structure) :: monolis
    integer(c_int), intent(in), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), value :: myrank, comm, commsize
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), target :: index(0:NP)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(0:recv_n_neib), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(0:send_n_neib), send_item(send_nitem)
    real(c_double), intent(in), target :: A(NDOF*NDOF*NZ)
    real(c_double), intent(in), target :: X(NDOF*NP)
    real(c_double), intent(in), target :: Y(NDOF*NP)

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NZ = NZ
    monolis%MAT%NDOF = NDOF
    monolis%MAT%A => A
    monolis%MAT%X => X
    monolis%MAT%B => Y
    monolis%MAT%index => index
    monolis%MAT%item => item

    !> for monoCOM
    monoliS%COM%myrank = myrank
    monoliS%COM%comm = comm
    monoliS%COM%commsize = commsize
    monoliS%COM%recv_n_neib = recv_n_neib
    monoliS%COM%recv_neib_pe => recv_neib_pe
    monoliS%COM%recv_index => recv_index
    monoliS%COM%recv_item => recv_item
    monoliS%COM%send_n_neib = send_n_neib
    monoliS%COM%send_neib_pe => send_neib_pe
    monoliS%COM%send_index => send_index
    monoliS%COM%send_item => send_item

    call monolis_matvec_product(monolis, monolis%MAT%X, monolis%MAT%B)
  end subroutine monolis_matvec_product_c

  subroutine monolis_inner_product_c(N, NDOF, X, Y, sum, comm) &
    & bind(c, name = "monolis_inner_product_c_main")
    use mod_monolis_com
    use mod_monolis_linalg
    implicit none
    integer(c_int), intent(in), value :: N, NDOF, comm
    real(c_double), target :: X(NDOF*N), Y(NDOF*N)
    real(c_double), value :: sum
    real(c_double) :: tspmv, tcomm
    type(monolis_com) :: monoCOM

    monoCOM%comm = comm
    call monolis_inner_product_R(monoCOM, n, ndof, X, Y, sum, tspmv, tcomm)
  end subroutine monolis_inner_product_c

  function monolis_allreduce_double_scalar_c(val, tag, comm) &
    & bind(c, name = "monolis_allreduce_double_scalar_c_main")
    use mod_monolis_linalg_com
    implicit none
    real(c_double) :: monolis_allreduce_double_scalar_c
    integer(c_int), intent(in), value :: tag, comm
    real(c_double), value :: val
    call  monolis_allreduce_R1(val, tag, comm)
    monolis_allreduce_double_scalar_c = val
  end function monolis_allreduce_double_scalar_c

  !> std lib
  subroutine monolis_qsort_int_c(array, iS, iE) &
    & bind(c, name = "monolis_qsort_int")
    implicit none
    integer(c_int), target :: array(1:iE-iS+1)
    integer(c_int), value :: iS, iE

    iS = iS
    iE = iE
    call monolis_qsort_int(array, iS, iE)
  end subroutine monolis_qsort_int_c

end module mod_monolis_wrapper
