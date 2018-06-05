module mod_monolis_linalg_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg_com
  implicit none

contains

  subroutine monolis_update_R(monoCOM, ndof, X, tcomm)
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_linalg_com
    implicit none
    type(monolis_com) :: monoCOM
    integer(kind=kint) :: i, ndof, ns, nr
    real(kind=kdouble) :: X(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm
    real(kind=kdouble), allocatable :: ws(:), wr(:)

    if( monoCOM%n_neib == 0 ) return

    ns = monoCOM%send_index(monoCOM%n_neib)
    nr = monoCOM%recv_index(monoCOM%n_neib)

    allocate(ws(ndof*ns), wr(ndof*nr))

    call monolis_SendRecv_R(monoCOM%n_neib, monoCOM%neib_pe, &
      &   monoCOM%send_index, monoCOM%send_item, &
      &   monoCOM%recv_index, monoCOM%recv_item, &
      &   ws, wr, X, ndof, monoCOM%comm)

    deallocate(ws, wr)
  end subroutine monolis_update_R

  subroutine monolis_update_I(monoCOM, ndof, X, tcomm)
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_linalg_com
    implicit none
    type(monolis_com) :: monoCOM
    integer(kind=kint) :: i, ndof, ns, nr
    integer(kind=kint) :: X(:)
    integer(kind=kint), allocatable :: ws(:), wr(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

    if( monoCOM%n_neib == 0 ) return

    ns = monoCOM%send_index(monoCOM%n_neib)
    nr = monoCOM%recv_index(monoCOM%n_neib)

    allocate(ws(ndof*ns), wr(ndof*nr))

    call monolis_SendRecv_I(monoCOM%n_neib, monoCOM%neib_pe, &
      &   monoCOM%send_index, monoCOM%send_item, &
      &   monoCOM%recv_index, monoCOM%recv_item, &
      &   ws, wr, X, ndof, monoCOM%comm)

    deallocate(ws, wr)
  end subroutine monolis_update_I

end module mod_monolis_linalg_util
