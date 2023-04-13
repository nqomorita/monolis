module mod_monolis_graph
  use mod_monolis_prm
  use mod_monolis_mesh
  use mod_monolis_util
  use mod_monolis_util_debug

  implicit none

contains

  subroutine monolis_part_graph(mesh, graph, n_domain)
    use iso_c_binding
    implicit none
    type(monolis_mesh) :: mesh
    type(monolis_graph) :: graph
    integer(kint), pointer :: part_id(:)
    integer(c_int), pointer :: index(:), item(:)
    integer(kint) :: i, n_domain

    call monolis_debug_header("monolis_part_graph")

    allocate(graph%node_domid_raw(mesh%nnode), source = 1)
    allocate(part_id(mesh%nnode), source = 0)

    if(n_domain == 1) return

    call monolis_convert_mesh_to_connectivity &
     & (mesh%nelem, mesh%nbase_func, mesh%elem, graph%ebase_func, graph%connectivity)

    call monolis_convert_connectivity_to_nodal_graph &
     & (mesh%nnode, mesh%nelem, graph%ebase_func, graph%connectivity, index, item)

    call monolis_get_partitioned_graph(mesh%nnode, index, item, n_domain, part_id)

    graph%nnode = mesh%nnode
    graph%nelem = mesh%nelem
    do i = 1, mesh%nnode
      graph%node_domid_raw(i) = part_id(i) + 1
    enddo

    deallocate(part_id)
  end subroutine monolis_part_graph

  subroutine monolis_part_nodal_graph(graph, n_domain)
    use iso_c_binding
    implicit none
    type(monolis_graph) :: graph
    integer(kint), pointer :: part_id(:)
    integer(c_int), pointer :: index(:), item(:)
    integer(kint) :: i, n_domain

    call monolis_debug_header("monolis_part_nodal_graph")

    allocate(graph%node_domid_raw(graph%N), source = 1)
    allocate(part_id(graph%N), source = 0)

    if(n_domain == 1) return

    call monolis_get_partitioned_graph(graph%N, graph%index, graph%item, n_domain, part_id)

    do i = 1, graph%N
      graph%node_domid_raw(i) = part_id(i) + 1
    enddo

    deallocate(part_id)
  end subroutine monolis_part_nodal_graph

  subroutine monolis_get_graph_from_graph_format(format, graph)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_graph_format) :: format
    integer(kint) :: i, j, in, nz

    graph%N = format%n_point
    graph%nnode = format%n_point

    nz = 0
    do i = 1, format%n_point
      nz = nz + format%n_adjacent(i)
    enddo

    allocate(graph%index(format%n_point+1), source = 0)
    do i = 1, format%n_point
      graph%index(i+1) = graph%index(i) + format%n_adjacent(i)
    enddo

    allocate(graph%item(nz), source = 0)
    graph%item = format%adjacent_id

    if(minval(format%point_id) == 0) graph%item = graph%item + 1
  end subroutine monolis_get_graph_from_graph_format

  subroutine monolis_convert_mesh_to_connectivity(nelem, nb, elem, ebase_func, connectivity)
    implicit none
    integer(kint) :: i, j, nb
    integer(kint) :: nelem, elem(:,:)
    integer(kint), pointer :: ebase_func(:), connectivity(:)

    allocate(ebase_func(nelem+1), source = 0)
    allocate(connectivity(nelem*nb), source = 0)

    do i = 1, nelem
      ebase_func(i+1) = i*nb
    enddo

    do i = 1, nelem
      do j = 1, nb
        connectivity(nb*(i-1) + j) = elem(j,i)
      enddo
    enddo
  end subroutine monolis_convert_mesh_to_connectivity

  subroutine monolis_convert_connectivity_to_nodal_graph &
    & (nnode, nelem, ebase_func, connectivity, index, item)
    use iso_c_binding
    implicit none
    integer(kint) :: nnode, nelem
    integer(kint), pointer :: ebase_func(:), connectivity(:)
    integer(c_int), pointer :: index(:), item(:)

    call monolis_convert_connectivity_to_nodal &
     & (nnode, nelem, ebase_func, connectivity, index, item)
  end subroutine monolis_convert_connectivity_to_nodal_graph

  subroutine monolis_convert_connectivity_to_nodal &
    & (nnode, nelem, ebase_func, connectivity, index, item)
    use iso_c_binding
    implicit none
    integer(kint) :: i, j, nnode, numflag
    integer(kint) :: nelem
    integer(kint), pointer :: ebase_func(:), connectivity(:)
    integer(c_int), pointer :: index(:), item(:)
    type(c_ptr) :: xadj, adjncy
#if WITH_METIS64
    integer(c_int64_t) :: nelem8, nnode8, numflag8
    integer(c_int64_t), pointer :: ebase_func8(:), connectivity8(:)
    integer(c_int64_t), pointer :: index8(:), item8(:)
#endif

    call monolis_debug_header("monolis_convert_connectivity_to_nodal")

    numflag = 0

    !> convert to 0 origin
    connectivity = connectivity - 1

#ifdef WITH_METIS
    call METIS_MESHTONODAL(nelem, nnode, ebase_func, connectivity, numflag, xadj, adjncy)
    call c_f_pointer(xadj, index, shape=[nnode+1])
    call c_f_pointer(adjncy, item, shape=[index(nnode+1)])
#elif WITH_METIS64
    nnode8 = nnode
    nelem8 = nelem
    numflag8 = numflag
    allocate(ebase_func8(nelem+1))
    allocate(connectivity8(ebase_func(nelem+1)))
    ebase_func8 = ebase_func
    connectivity8 = connectivity
    call METIS_MESHTONODAL(nelem8, nnode8, ebase_func8, connectivity8, numflag8, xadj, adjncy)
    call c_f_pointer(xadj, index8, shape=[nnode+1])
    call c_f_pointer(adjncy, item8, shape=[index8(nnode+1)])
    allocate(index(nnode+1))
    allocate(item(index8(nnode+1)))
    index = index8
    item = item8
#else
    call monolis_warning_header("monolis_convert_connectivity_to_nodal: METIS is NOT enabled")
    stop
#endif

    !> convert to 1 origin
    connectivity = connectivity + 1

    !> convert to 1 origin
    item = item + 1
  end subroutine monolis_convert_connectivity_to_nodal

  subroutine monolis_get_partitioned_graph(nnode, index, item, n_domain, part_id)
    use iso_c_binding
    implicit none
    integer(kint), pointer :: part_id(:)
    integer(c_int), pointer :: index(:), item(:)
    integer(c_int), pointer :: node_wgt(:) => null()
    integer(c_int), pointer :: edge_wgt(:) => null()
    integer(kint) :: n_domain, nnode

    call monolis_debug_header("monolis_get_partitioned_graph")
    call monolis_get_mesh_part_kway(nnode, index, item, n_domain, node_wgt, edge_wgt, part_id)
  end subroutine monolis_get_partitioned_graph

  subroutine monolis_get_partitioned_graph_with_node_weight(nnode, index, item, n_domain, wgt, part_id)
    use iso_c_binding
    implicit none
    integer(kint), pointer :: part_id(:), wgt(:)
    integer(c_int), pointer :: index(:), item(:)
    integer(c_int), pointer :: node_wgt(:)
    integer(c_int), pointer :: edge_wgt(:) => null()
    integer(kint) :: n_domain, nnode

    call monolis_debug_header("monolis_get_partitioned_graph_with_node_weight")
    allocate(node_wgt(nnode), source = 0)
    node_wgt = wgt

    call monolis_get_mesh_part_kway(nnode, index, item, n_domain, node_wgt, edge_wgt, part_id)

    deallocate(node_wgt)
  end subroutine monolis_get_partitioned_graph_with_node_weight

  subroutine monolis_get_mesh_part_kway(nnode, index, item, npart, node_wgt, edge_wgt, part_id)
    use iso_c_binding
    implicit none
    integer(kint) :: nnode, ncon, npart, objval
    integer(kint), pointer :: part_id(:)
    integer(c_int), pointer :: node_wgt(:)
    integer(c_int), pointer :: edge_wgt(:)
    integer(kint), pointer :: vsize(:) => null()
    integer(kint), pointer :: ubvec(:) => null()
    real(kdouble), pointer :: options(:) => null()
    real(kdouble), pointer :: tpwgts(:) => null()
    integer(c_int), pointer :: index(:), item(:)
#if WITH_METIS64
    integer(c_int64_t) :: nnode8, ncon8, npart8, objval8
    integer(c_int64_t), pointer :: part_id8(:)
    integer(c_int64_t), pointer :: node_wgt8(:)
    integer(c_int64_t), pointer :: edge_wgt8(:) => null()
    integer(c_int64_t), pointer :: vsize8(:) => null()
    integer(c_int64_t), pointer :: ubvec8(:)  => null()
    integer(c_int64_t), pointer :: index8(:), item8(:)
#endif

    call monolis_debug_header("monolis_get_mesh_part_kway")

    if(npart /= 1)then
      ncon = 1
      !> convert to 0 origin
      item = item - 1

#ifdef WITH_METIS
      call METIS_PARTGRAPHRECURSIVE(nnode, ncon, index, item, node_wgt, vsize, edge_wgt, npart, tpwgts, ubvec, &
      !call METIS_PARTGRAPHKWAY(nnode, ncon, index, item, node_wgt, vsize, edge_wgt, npart, tpwgts, ubvec, &
        & options, objval, part_id)
#elif WITH_METIS64
      ncon8 = 1
      npart8 = npart
      nnode8 = nnode
      allocate(node_wgt8(nnode))
      node_wgt8 = node_wgt
      allocate(index8(nnode+1))
      index8 = index
      allocate(item8(index(nnode+1)))
      item8 = item
      allocate(part_id8(nnode))
      call METIS_PARTGRAPHRECURSIVE(nnode8, ncon8, index8, item8, node_wgt8, vsize8, edge_wgt8, npart8, tpwgts, ubvec8, &
      !call METIS_PARTGRAPHKWAY(nnode, ncon, index, item, node_wgt, vsize, edge_wgt, npart, tpwgts, ubvec, &
        & options, objval8, part_id8)
      part_id = part_id8
#else
    call monolis_warning_header("monolis_get_mesh_part_kway: METIS is NOT enabled")
    stop
#endif

      !> convert to 1 origin
      item = item + 1
    endif
  end subroutine monolis_get_mesh_part_kway

end module mod_monolis_graph
