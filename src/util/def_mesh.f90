module mod_monolis_mesh
  use mod_monolis_prm
  implicit none

  type monolis_mesh
    integer(kint) :: nnode, nnode_in
    integer(kint) :: nelem, nelem_in
    integer(kint) :: nbase_func
    integer(kint), allocatable :: nid(:)
    real(kdouble), allocatable :: node(:,:)
    integer(kint), allocatable :: eid(:)
    integer(kint), allocatable :: elem(:,:)
  end type monolis_mesh

  type monolis_graph
    integer(kint) :: N = 0
    integer(kint) :: nelem = 0
    !> node base
    integer(kint), pointer :: node_domid_raw(:) => null()
    integer(kint), pointer :: node_domid(:) => null()
    !> elem base
    integer(kint), pointer :: elem_domid(:) => null()
    integer(kint), pointer :: elem_domid_uniq(:) => null()
    !> elem base
    integer(kint), pointer :: ebase_func(:) => null()
    integer(kint), pointer :: connectivity(:) => null()
  end type monolis_graph

  type monolis_node_list
    integer(kint) :: nnode, nnode_in, nnode_out
    integer(kint) :: nelem, nelem_in, nelem_out
    integer(kint), allocatable :: nid(:)
    integer(kint), allocatable :: nid_perm(:)
    integer(kint), allocatable :: eid(:)
    integer(kint), allocatable :: recv_item_perm(:)
    integer(kint), allocatable :: recv_item_domid(:)
  end type monolis_node_list

contains

  subroutine monolis_get_connectivity(mesh, eid, nbase_func, conn)
    implicit none
    type(monolis_mesh) :: mesh
    integer(kint) :: i, eid, nbase_func, conn(nbase_func)

    do i = 1, nbase_func
      conn(i) = mesh%elem(i,eid)
    enddo
  end subroutine monolis_get_connectivity
end module mod_monolis_mesh