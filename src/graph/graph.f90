module mod_monolis_graph
  use mod_monolis_prm
  use mod_monolis_util

  implicit none

contains

  subroutine monolis_get_mesh_to_nodal(nnode, nelem, nbase, elem, index, item)
    implicit none
    integer(kint) :: i, j, nnode, nbase, nind, numflag
    integer(kint) :: nelem, elem(:,:)
    integer(kint), pointer :: eptr(:) => null()
    integer(kint), pointer :: eind(:) => null()
    integer(kint), pointer :: index(:), item(:)

    nind = nelem*nbase
    numflag = 1
    allocate(eptr(nelem+1), source = 0)
    allocate(eind(nind), source = 0)

    do i = 1, nelem
      eptr(i+1) = i*nbase
    enddo

    do i = 1, nelem
      do j = 1, nbase
        eind(nbase*(i-1) + j) = elem(j,i) - 1
      enddo
    enddo

#ifdef WITH_METIS
    call METIS_MESHTONODAL(nelem, nnode, eptr, eind, numflag, index, item)
#endif
  end subroutine monolis_get_mesh_to_nodal

  subroutine monolis_get_mesh_part_kway(nnode, index, item, npart, part_id)
    implicit none
    integer(kint) :: nnode, ncon, npart, objval
    integer(kint), pointer :: part_id(:)
    integer(kint), pointer :: eptr(:)   => null()
    integer(kint), pointer :: eind(:)   => null()
    integer(kint), pointer :: vwgtm(:)  => null()
    integer(kint), pointer :: vsize(:)  => null()
    integer(kint), pointer :: adjwgt(:) => null()
    integer(kint), pointer :: ubvec(:)  => null()
    real(kdouble), pointer :: options(:) => null()
    real(kdouble), pointer :: tpwgts(:) => null()
    integer(kint), pointer :: index(:), item(:)

    allocate(part_id(nnode), source = 0)

    if(npart /= 1)then
      ncon = 1
      call METIS_PARTGRAPHKWAY(nnode, ncon, index, item, vwgtm, vsize, adjwgt, npart, tpwgts, ubvec, &
        & options, objval, part_id)
    endif
  end subroutine monolis_get_mesh_part_kway

end module mod_monolis_graph
