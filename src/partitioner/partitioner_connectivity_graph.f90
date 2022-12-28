program monolis_partitioner_nodal_graph
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_graph
  use mod_monolis_comm_overlap
  use mod_monolis_io
  use mod_monolis_io_arg
  use mod_monolis_comm_util
  implicit none
  type(monolis_graph) :: graph
  type(monolis_graph_format) :: graph_format
  type(monolis_mesh), allocatable :: mesh(:)
  type(monolis_node_list), allocatable :: node_list(:)
  integer(kint) :: n_domain, nid, shift
  character :: fname*100

  call monolis_global_initialize()

  call monolis_get_connectivity_part_arg(fname, n_domain)

  if(n_domain <= 1) stop

  call monolis_input_graph(fname, graph_format)

  call monolis_get_graph_from_graph_format(graph_format, graph)

  allocate(mesh(n_domain))

  call monolis_par_input_node_id(n_domain, mesh)

  call monolis_par_input_node_n_internal(n_domain, mesh)

  allocate(node_list(n_domain))

  call get_elem_domid_uniq(graph, n_domain, mesh, shift)

  do nid =  1, n_domain
    call get_connectivity_at_subdomain(graph, mesh(nid), node_list(nid), nid, shift)
  enddo

  call monolis_output_parted_connectivity_graph(fname, mesh, graph, graph_format, node_list, n_domain)

  call monolis_part_finalize()

  call monolis_global_finalize()

contains

  subroutine get_elem_domid_uniq(graph, n_domain, mesh, shift)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_mesh) :: mesh(:)
    integer(kint) :: n_domain
    integer(kint) :: i, j, in, nnode_all, nelem, jS, jE, shift

    nnode_all = 0
    shift = 0
    do i = 1, n_domain
      nnode_all = nnode_all + mesh(i)%nnode_in
      if(minval(mesh(i)%nid) == 0) shift = 1
    enddo

    allocate(graph%node_domid_raw(nnode_all), source = 0)

    nnode_all = 0
    do i = 1, n_domain
      do j = 1, mesh(i)%nnode_in
        in = mesh(i)%nid(j) + shift
        graph%node_domid_raw(in) = i
      enddo
    enddo

    nelem = graph%nnode
    allocate(graph%elem_domid_uniq(nelem), source = 2147483647)
    do i = 1, nelem
      jS = graph%index(i) + 1
      jE = graph%index(i+1)
      do j = jS, jE
        in = graph%item(j)
        if(graph%node_domid_raw(in) < graph%elem_domid_uniq(i))then
          graph%elem_domid_uniq(i) = graph%node_domid_raw(in)
        endif
      enddo
    enddo
  end subroutine get_elem_domid_uniq

  subroutine get_connectivity_at_subdomain(graph, mesh, local_node, nid, shift)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_mesh) :: mesh
    type(monolis_node_list) :: local_node
    integer(kint) :: nelem, nnode_all, i, j, jS, jE, jn, nid, count_in, count_out, in, shift
    logical, allocatable :: is_node_used(:), is_in(:)

    nnode_all = maxval(graph%item)
    allocate(is_node_used(nnode_all), source = .false.)

    do i = 1, mesh%nnode
      in = mesh%nid(i) + shift
      if(nnode_all < in) cycle
      is_node_used(in) = .true.
    enddo

    nelem = graph%nnode
    allocate(is_in(nelem), source = .false.)

    do i = 1, nelem
      jS = graph%index(i) + 1
      jE = graph%index(i+1)
      jn = 0
      do j = jS, jE
        in = graph%item(j)
        if(is_node_used(in)) jn = jn + 1
      enddo
      if(jn == jE-jS+1) is_in(i) = .true.
    enddo

    in = 0
    do j = 1, nelem
      if(is_in(j)) in = in + 1
    enddo

    local_node%nelem = in
    allocate(local_node%eid(in), source = 0)

    !> internal element (global unique element)
    in = 0
    count_in = 0
    do j = 1, nelem
      if(is_in(j) .and. graph%elem_domid_uniq(j) == nid)then
        in = in + 1
        count_in = count_in + 1
        local_node%eid(in) = j
      endif
    enddo

    !> external element
    count_out = 0
    do j = 1, nelem
      if(is_in(j) .and. graph%elem_domid_uniq(j) /= nid)then
        in = in + 1
        count_out = count_out + 1
        local_node%eid(in) = j
      endif
    enddo

    local_node%nelem_in = count_in
    local_node%nelem_out = count_out
  end subroutine get_connectivity_at_subdomain

  subroutine monolis_part_finalize()
    implicit none
    if(allocated(node_list)) deallocate(node_list)
  end subroutine monolis_part_finalize

end program monolis_partitioner_nodal_graph
