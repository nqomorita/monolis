module mod_monolis_linalg_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_linalg_com
  implicit none

contains

  subroutine monolis_update_R(monoCOM, ndof, X, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: ndof, ns, nr
    real(kdouble) :: X(:)
    real(kdouble) :: t1, t2, tcomm
    real(kdouble), allocatable :: ws(:), wr(:)

    if(monoCOM%send_n_neib == 0 .and. monoCOM%recv_n_neib == 0) return
    if(.not. monoCOM%is_overlap) return

    ns = monoCOM%send_index(monoCOM%send_n_neib)
    nr = monoCOM%recv_index(monoCOM%recv_n_neib)

    allocate(ws(ndof*ns), wr(ndof*nr))

    t1 = monolis_get_time()
    call monolis_SendRecv_pre_R(monoCOM%send_n_neib, monoCOM%send_neib_pe, &
       & monoCOM%recv_n_neib, monoCOM%recv_neib_pe, &
       & monoCOM%send_index, monoCOM%send_item, &
       & monoCOM%recv_index, monoCOM%recv_item, &
       & ws, wr, X, ndof, monoCOM%comm)
    t2 = monolis_get_time()
    tcomm = tcomm + t2 - t1

    deallocate(ws, wr)
  end subroutine monolis_update_R

  subroutine monolis_update_I(monoCOM, ndof, X, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: ndof, ns, nr
    integer(kint) :: X(:)
    integer(kint), allocatable :: ws(:), wr(:)
    real(kdouble), optional :: tcomm

    if(monoCOM%send_n_neib == 0 .and. monoCOM%recv_n_neib == 0) return

    ns = monoCOM%send_index(monoCOM%send_n_neib)
    nr = monoCOM%recv_index(monoCOM%recv_n_neib)

    allocate(ws(ndof*ns), wr(ndof*nr))

    call monolis_SendRecv_pre_I(monoCOM%send_n_neib, monoCOM%send_neib_pe, &
       & monoCOM%recv_n_neib, monoCOM%recv_neib_pe, &
       & monoCOM%send_index, monoCOM%send_item, &
       & monoCOM%recv_index, monoCOM%recv_item, &
       & ws, wr, X, ndof, monoCOM%comm)

    deallocate(ws, wr)
  end subroutine monolis_update_I
end module mod_monolis_linalg_util
