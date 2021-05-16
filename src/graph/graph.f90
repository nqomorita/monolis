module mod_monolis_graph
  use mod_monolis_prm
  use mod_monolis_mesh
  use mod_monolis_util

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

    call monolis_convert_connectivity_to_nodal &
     & (mesh%nnode, mesh%nelem, graph%ebase_func, graph%connectivity, index, item)

    call monolis_get_mesh_part_kway(mesh%nnode, index, item, n_domain, part_id)

    graph%nelem = mesh%nelem
    do i = 1, mesh%nnode
      graph%node_domid_raw(i) = part_id(i) + 1
    enddo

    deallocate(part_id)
  end subroutine monolis_part_graph

  subroutine monolis_convert_mesh_to_connectivity(nelem, nb, elem, ebase_func, connectivity)
    use iso_c_binding
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

  subroutine monolis_convert_connectivity_to_nodal &
    & (nnode, nelem, ebase_func, connectivity, index, item)
    use iso_c_binding
    implicit none
    integer(kint) :: i, j, nnode, numflag
    integer(kint) :: nelem
    integer(kint), pointer :: ebase_func(:), connectivity(:)
    integer(c_int), pointer :: index(:), item(:)
    type(c_ptr) :: xadj, adjncy

    call monolis_debug_header("monolis_convert_connectivity_to_nodal")

    numflag = 0

    !> convert to 0 origin
    connectivity = connectivity - 1

#ifdef WITH_METIS
    call METIS_MESHTONODAL(nelem, nnode, ebase_func, connectivity, numflag, xadj, adjncy)
    call c_f_pointer(xadj, index, shape=[nnode+1])
    call c_f_pointer(adjncy, item, shape=[index(nnode+1)])
#else
    call monolis_warning_header("monolis_convert_connectivity_to_nodal: METIS is NOT enabled")
    stop
#endif

    !> convert to 1 origin
    item = item + 1
  end subroutine monolis_convert_connectivity_to_nodal

  subroutine monolis_get_mesh_part_kway(nnode, index, item, npart, part_id)
    use iso_c_binding
    implicit none
    integer(kint) :: nnode, ncon, npart, objval
    integer(kint), pointer :: part_id(:)
    integer(kint), pointer :: vwgtm(:)  => null()
    integer(kint), pointer :: vsize(:)  => null()
    integer(kint), pointer :: adjwgt(:) => null()
    integer(kint), pointer :: ubvec(:)  => null()
    real(kdouble), pointer :: options(:) => null()
    real(kdouble), pointer :: tpwgts(:) => null()
    integer(c_int), pointer :: index(:), item(:)

    call monolis_debug_header("monolis_get_mesh_part_kway")

    if(npart /= 1)then
      ncon = 1
      !> convert to 0 origin
      item = item - 1
#ifdef WITH_METIS
      call METIS_PARTGRAPHRECURSIVE(nnode, ncon, index, item, vwgtm, vsize, adjwgt, npart, tpwgts, ubvec, &
      !call METIS_PARTGRAPHKWAY(nnode, ncon, index, item, vwgtm, vsize, adjwgt, npart, tpwgts, ubvec, &
        & options, objval, part_id)
      !> convert to 1 origin
      item = item + 1
#else
    call monolis_warning_header("monolis_get_mesh_part_kway: METIS is NOT enabled")
    stop
#endif
    endif
  end subroutine monolis_get_mesh_part_kway

end module mod_monolis_graph
