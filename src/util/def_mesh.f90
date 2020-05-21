module mod_monolis_mesh
  use mod_monolis_prm
  implicit none

  type monolis_mesh
    integer(kint) :: n_node
    integer(kint) :: nnode_in, nnode_out
    integer(kint), allocatable :: nid(:)
    real(kdouble), allocatable :: node(:,:)
    integer(kint) :: n_elem
    integer(kint), allocatable :: elem(:,:)
    integer(kint) :: n_base_func
  end type monolis_mesh

  type monolis_node_list
    integer(kint) :: n_node
    integer(kint) :: n_node_in, n_node_out
    integer(kint) :: n_elem
    integer(kint), allocatable :: nid(:)
    integer(kint), allocatable :: nid_perm(:)
    integer(kint), allocatable :: eid(:)
  end type monolis_node_list

end module mod_monolis_mesh