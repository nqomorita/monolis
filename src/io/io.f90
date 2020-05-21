module mod_monolis_io
  use mod_monolis_prm
  use mod_monolis_util
  implicit none

contains

  subroutine monolis_input_mesh_node(fname, nnode, node, nid)
    implicit none
    integer(kint) :: nnode, i
    real(kdouble), allocatable :: node(:,:)
    integer(kint), optional, allocatable :: nid(:)
    character :: fname*100

    open(20, file = fname, status = "old")
      read(20,*) nnode
      call monolis_debug_int("nnode", nnode)

      allocate(node(3,nnode), source = 0.0d0)

      if(present(nid))then
        allocate(nid(nnode), source = 0)
        do i = 1, nnode
          read(20,*) nid(i), node(1,i), node(2,i), node(3,i)
        enddo
      else
        do i = 1, nnode
          read(20,*) node(1,i), node(2,i), node(3,i)
        enddo
      endif
    close(20)
  end subroutine monolis_input_mesh_node

  subroutine monolis_input_mesh_elem(fname, nelem, nbase, elem, eid)
    implicit none
    integer(kint) :: nelem, nbase, i, j
    integer(kint), allocatable :: elem(:,:)
    integer(kint), optional, allocatable :: eid(:)
    character :: fname*100

    open(20, file = fname, status = "old")
      read(20,*) nelem, nbase
      call monolis_debug_int("nelem", nelem)
      call monolis_debug_int("nbase_func", nbase)

      allocate(elem(nbase,nelem), source = 0)

      if(present(eid))then
        allocate(eid(nelem), source = 0)
        do i = 1, nelem
          read(20,*) eid(i), (elem(j,i), j = 1, nbase)
        enddo
      else
        do i = 1, nelem
          read(20,*) (elem(j,i), j = 1, nbase)
        enddo
      endif
    close(20)
  end subroutine monolis_input_mesh_elem

  subroutine monolis_output_mesh_node(fname, nnode, nnode_in, nid, node)
    implicit none
    integer(kint) :: nnode, nnode_in
    integer(kint) :: i, in, nid(:)
    real(kdouble) :: node(:,:)
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0,x,i0)")nnode_in, nnode
      do i = 1, nnode
        in = nid(i)
        write(20,"(1p3e22.14)") node(1,in), node(2,in), node(3,in)
      enddo
    close(20)
  end subroutine monolis_output_mesh_node

  subroutine monolis_output_mesh_elem(fname, nelem, nbase, eid, elem, perm)
    implicit none
    integer(kint) :: i, j, in, jn, nelem, nbase, eid(:), elem(:,:), perm(:)
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0,x,i0)")nelem, nbase
      do i = 1, nelem
        jn = eid(i)
        do j = 1, nbase
          in = elem(j,jn)
          write(20,"(x,i0,$)") perm(in)
        enddo
        write(20,*)""
      enddo
    close(20)
  end subroutine monolis_output_mesh_elem

  subroutine monolis_output_mesh_global_nid(fname, nnode, global_nid, nid)
    implicit none
    integer(kint) :: i, in, nnode, global_nid(:), nid(:)
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0)")nnode
      do i = 1, nnode
        in = nid(i)
        write(20,"(i0)") global_nid(in)
      enddo
    close(20)
  end subroutine monolis_output_mesh_global_nid

  subroutine monolis_output_mesh_global_eid(fname, nelem, global_eid, eid)
    implicit none
    integer(kint) :: i, in, nelem, global_eid(:), eid(:)
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0,x,i0)")nelem
      do i = 1, nelem
        in = eid(i)
        write(20,"(i0,x,i0)") global_eid(in)
      enddo
    close(20)
  end subroutine monolis_output_mesh_global_eid

  subroutine monolis_output_mesh_comm(fname, n_neib, neib_pe, index, item)
    implicit none
    integer(kint) :: n_neib
    integer(kint) :: neib_pe(:)
    integer(kint) :: index(0:n_neib)
    integer(kint) :: item(:)
    integer(kint) :: i
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0,x,i0)")n_neib, index(n_neib)
      do i = 1, n_neib
        write(20,"(i0)")neib_pe(i) - 1 !> for MPI
      enddo

      do i = 0, n_neib
        write(20,"(i0)")index(i)
      enddo
      do i = 1, index(n_neib)
        write(20,"(i0)")item(i)
      enddo

      write(20,"(i0,x,i0)")0, 0 !> for non-over
    close(20)
  end subroutine monolis_output_mesh_comm

  subroutine monolis_visual_parted_mesh(nnode, node, nelem, nbase, elem)
    implicit none
    integer(kint) :: i, j, in
    integer(kint) :: nnode, nelem, nbase, elem(:,:)
    real(kdouble) :: node(:,:)
    character :: etype*6, output_dir*100, fname*100

    call monolis_debug_header("monolis_visual_parted_mesh")

    output_dir = "visual/"
    call system('if [ ! -d visual ]; then (echo "** create visual"; mkdir -p visual); fi')

    open(20, file = fname, status = "replace")
      write(20,"(5i12)") nnode, nelem, 0, 0, 0
      do i = 1, nnode
        write(20,"(i0,1p3e12.5)") i, node(1,i), node(2,i), node(3,i)
      enddo

      if(nbase == 4) etype = " tet  "
      if(nbase == 8) etype = " hex  "
      if(nbase ==10) etype = " tet2 "

      do i = 1, nelem
        write(20,"(i0,i4,a,$)") i, 0, etype
        do j = 1, nbase
          write(20,"(i12,$)") elem(j,i)
        enddo
        write(20,*)""
      enddo

      write(20,"(a)")"1 1"
      write(20,"(a)")"node_domid, unknown"
      do i = 1, nnode
        !write(20,"(i0,x,i0,x,i0,x,i0)") i, bp_graph%node_domid_raw(i)!, bp_graph%elem_domid(i), in
      enddo

      write(20,"(a)")"1 1"
      write(20,"(a)")"elem_domid, unknown"
      do i = 1, nelem
        !write(20,"(i0,x,i0,x,i0,x,i0)") i, bp_graph%elem_domid_raw(i), bp_graph%elem_domid(i), in
      enddo
    close(20)
  end subroutine monolis_visual_parted_mesh

end module mod_monolis_io
