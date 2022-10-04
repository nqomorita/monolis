module mod_monolis_graph_comm
  use mod_monolis_prm
  use mod_monolis_mesh
  use mod_monolis_util

  implicit none
  private
  public :: monolis_get_comm_table

!  type id_list
!    integer(kint) :: n
!    integer(kint), allocatable :: id(:)
!  end type id_list

contains

  subroutine monolis_get_comm_table(monolis, n_id, global_node_id)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_id, global_node_id(:)
    integer(kint) :: commsize
    type(id_list), allocatable :: nids(:), recv_array(:)

!    if(monolis%COM%commsize <= 1)then
!      (monolis%COM%commsize = 1
!      return
!    endif
!
!    if(n_id /= monolis%MAT%N)then
!      stop "monolis_get_comm_table"
!    endif
!
!    !> set global node id
!    allocate(monolis%COM%global_node_id(n_id), source = 0)
!    monolis%COM%global_node_id = global_node_id
!
!    commsize = monolis%COM%commsize
!    allocate(nids(commsize))
!    allocate(recv_array(commsize))
!    call monolis_allgather_I1(n_id, recv_array, monolis%COM%comm)
!
!    do i = 1, commsize
!      if()then
!
!      endif
!    enddo
!
!!monolis_Isend_I(n, ws, pe_id, comm, req)
!!monolis_Irecv_I(n, ws, pe_id, comm, req)
!
!    !integer(kint)          :: internal_nnode
!
!    !integer(kint)          :: recv_n_neib
!    !integer(kint), pointer :: recv_neib_pe(:) => null()
!    !integer(kint), pointer :: recv_index(:)   => null()
!    !integer(kint), pointer :: recv_item(:)    => null()
!
!    !integer(kint)          :: send_n_neib
!    !integer(kint), pointer :: send_neib_pe(:) => null()
!    !integer(kint), pointer :: send_index(:)   => null()
!    !integer(kint), pointer :: send_item(:)    => null()
  end subroutine monolis_get_comm_table

end module mod_monolis_graph_comm
