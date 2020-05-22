module mod_monolis_comm_util
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  implicit none

contains

  subroutine get_overlap_domain(mesh, graph, n_domain)
    implicit none
    type(monolis_mesh) :: mesh
    type(monolis_graph) :: graph
    integer(kint) :: nelem, n_domain
    integer(kint) :: i, j, in, id, id1, maxid
    logical :: is_overlap

    call monolis_debug_header("get_overlap_domain")

    nelem = mesh%nelem
    allocate(graph%elem_domid(nelem), source = 0)
    !allocate(graph%elem_domid_raw(nelem), source = 0)
    allocate(graph%elem_domid_uniq(nelem), source = 0)

    if(n_domain == 1)then
      graph%elem_domid = 1
      graph%elem_domid_raw = 1
      graph%elem_domid_uniq = 1
      return
    endif

    do i = 1, nelem
      in = mesh%elem(1,i)
      id1 = graph%node_domid_raw(in)
      maxid = id1
      is_overlap = .false.
      do j = 2, mesh%nbase_func
        in = mesh%elem(j,i)
        id = graph%node_domid_raw(in)
        if(id /= id1) is_overlap = .true.
        maxid = max(id1, id)
      enddo
      graph%elem_domid_uniq(i) = maxid
      if(is_overlap)then
        graph%elem_domid(i) = -1
      else
        graph%elem_domid(i) = id1
      endif
    enddo

    !graph%elem_domid_raw = graph%elem_domid

    do i = 1, nelem
      if(graph%elem_domid(i) == 0) stop "get_overlap_domain"
    enddo
  end subroutine get_overlap_domain

  subroutine get_nnode_and_nid_at_subdomain(mesh, graph, local_node, nid, avg)
    implicit none
    type(monolis_mesh) :: mesh
    type(monolis_graph) :: graph
    type(monolis_node_list) :: local_node
    integer(kint) :: nnode, nid, in, j, jn, k
    integer(kint) :: count_in, count_out, avg
    logical, allocatable :: is_in(:)

    nnode = mesh%nnode
    allocate(is_in(nnode), source = .false.)

    do j = 1, local_node%nelem
      in = local_node%eid(j)
      do k = 1, mesh%nbase_func
        jn = mesh%elem(k,in)
        is_in(jn) = .true.
      enddo
    enddo

    count_in = 0
    count_out = 0
    do j = 1, nnode
      if(is_in(j))then
        if(graph%node_domid_raw(j) /= nid)then
          count_out = count_out + 1
        else
          count_in = count_in + 1
        endif
      endif
    enddo

    !> save to structure
    local_node%nnode     = count_in + count_out
    local_node%nnode_in  = count_in
    local_node%nnode_out = count_out

write(*,"(4i10)") nid, local_node%nnode, local_node%nnode_in, local_node%nnode_out

    avg = avg + local_node%nnode
    allocate(local_node%nid(local_node%nnode), source = 0)

    in = 1
    do j = 1, nnode
      if(is_in(j) .and. graph%node_domid_raw(j) == nid)then
        local_node%nid(in) = j
        in = in + 1
      endif
    enddo
    do j = 1, nnode
      if(is_in(j) .and. graph%node_domid_raw(j) /= nid)then
        local_node%nid(in) = j
        in = in + 1
      endif
    enddo
  end subroutine get_nnode_and_nid_at_subdomain

  subroutine get_nelem_and_eid_at_subdomain(mesh, graph, local_node, nid)
    implicit none
    type(monolis_mesh) :: mesh
    type(monolis_graph) :: graph
    type(monolis_node_list) :: local_node
    integer(kint) :: nelem, nid, in, j, k
    logical, allocatable :: is_in(:)

    call monolis_debug_header("get_nelem_and_eid_at_subdomain")

    nelem = mesh%nelem
    allocate(is_in(nelem), source = .false.)

    is_in = .false.
    do j = 1, nelem
      do k = 1, mesh%nbase_func
        in = mesh%elem(k,j)
        if(graph%node_domid_raw(in) == nid)then
          is_in(j) = .true.
        endif
      enddo
    enddo

    in = 0
    do j = 1, nelem
      if(is_in(j)) in = in + 1
    enddo

    local_node%nelem = in
    allocate(local_node%eid(in))

    in = 0
    do j = 1, nelem
      if(is_in(j))then
        in = in + 1
        local_node%eid(in) = j
      endif
    enddo
  end subroutine get_nelem_and_eid_at_subdomain

end module mod_monolis_comm_util
