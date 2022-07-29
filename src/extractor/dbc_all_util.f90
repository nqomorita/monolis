module mod_monolis_dbc_all_util
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  use mod_monolis_hash
  use mod_monolis_stdlib
  implicit none

  private
  public :: monolis_get_surf_node
  public :: output_dbc
  public :: output_surf

  integer(kint) :: i341(3,4), i361(4,6)
  type(type_monolis_hash_tree) :: hash_tree

contains

  subroutine output_dbc(foname, mesh, is_surf_node, n_block, val)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh
    integer(kint) :: i, j, in, shift, n_block
    integer(kint) :: is_surf_node(:)
    real(kdouble) :: val(:)
    character :: foname*100

    shift = 0
    if(minval(mesh%nid) == 0) shift = -1

    in = 0
    do i = 1, mesh%nnode
      if(is_surf_node(i) == 1) in = in + 1
    enddo

    open(20, file = trim(foname), status = "replace")
      write(20,*) n_block*in, n_block
      do i = 1, mesh%nnode
        do j = 1, n_block
          if(is_surf_node(i) == 1) write(20,*) i+shift, j+shift, val(j)
        enddo
      enddo
    close(20)
  end subroutine output_dbc

  subroutine output_surf(foname, mesh, nbase_func, nsurf, nsurf_node, is_open_surf)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh
    integer(kint) :: i, j, in, shift, conn(nbase_func), nbase_func, nsurf, nsurf_node
    integer(kint) :: is_open_surf(:,:), eid
    character :: foname*100

    i341(1,1) = 1; i341(2,1) = 2; i341(3,1) = 3
    i341(1,2) = 1; i341(2,2) = 2; i341(3,2) = 4
    i341(1,3) = 2; i341(2,3) = 3; i341(3,3) = 4
    i341(1,4) = 3; i341(2,4) = 1; i341(3,4) = 4

    i361(1,1) = 1; i361(2,1) = 2; i361(3,1) = 3; i361(4,1) = 4
    i361(1,2) = 5; i361(2,2) = 6; i361(3,2) = 7; i361(4,2) = 8
    i361(1,3) = 1; i361(2,3) = 2; i361(3,3) = 5; i361(4,3) = 6
    i361(1,4) = 2; i361(2,4) = 3; i361(3,4) = 6; i361(4,4) = 7
    i361(1,5) = 3; i361(2,5) = 4; i361(3,5) = 7; i361(4,5) = 8
    i361(1,6) = 4; i361(2,6) = 1; i361(3,6) = 8; i361(4,6) = 5

    shift = 0
    if(minval(mesh%nid) == 0) shift = -1

    in = 0
    do i = 1, mesh%nelem
      do j = 1, nbase_func
        if(is_open_surf(j,i) == 1) in = in + 1
      enddo
    enddo

    open(20, file = trim(foname), status = "replace")
      write(20,*) in

      do eid = 1, mesh%nelem
        call monolis_get_connectivity(mesh, eid, nbase_func, conn)
        do i = 1, nsurf
          if(is_open_surf(j,i) == 1)then
            do j = 1, nsurf_node
              if(nbase_func == 4)then
                in = conn(i341(j,i))
              elseif(nbase_func == 8)then
                in = conn(i361(j,i))
              endif
              write(20,*) in-shift
            enddo
            write(20,*)""
          endif
        enddo
      enddo
    close(20)
  end subroutine output_surf

  subroutine monolis_get_surf_node(mesh, nbase_func, nsurf, nsurf_node, is_surf_node, is_open_surf)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh
    integer(kint) :: nbase_func, conn(nbase_func), nsurf, nsurf_node
    integer(kint) :: i, j, in, eid, tmp
    character :: ckey*27
    logical :: is_exist, is_pushed
    integer(kint), allocatable :: is_inner(:,:), is_surf_node(:)
    integer(kint), optional, allocatable :: is_open_surf(:,:)

    i341(1,1) = 1; i341(2,1) = 2; i341(3,1) = 3
    i341(1,2) = 1; i341(2,2) = 2; i341(3,2) = 4
    i341(1,3) = 2; i341(2,3) = 3; i341(3,3) = 4
    i341(1,4) = 3; i341(2,4) = 1; i341(3,4) = 4

    i361(1,1) = 1; i361(2,1) = 2; i361(3,1) = 3; i361(4,1) = 4
    i361(1,2) = 5; i361(2,2) = 6; i361(3,2) = 7; i361(4,2) = 8
    i361(1,3) = 1; i361(2,3) = 2; i361(3,3) = 5; i361(4,3) = 6
    i361(1,4) = 2; i361(2,4) = 3; i361(3,4) = 6; i361(4,4) = 7
    i361(1,5) = 3; i361(2,5) = 4; i361(3,5) = 7; i361(4,5) = 8
    i361(1,6) = 4; i361(2,6) = 1; i361(3,6) = 8; i361(4,6) = 5

    call monolis_debug_int("nbase_func", nbase_func)
    call monolis_debug_int("nsurf", nsurf)
    call monolis_debug_int("nsurf_node", nsurf_node)

    call monolis_hash_init(hash_tree)

    do eid = 1, mesh%nelem
      call monolis_get_connectivity(mesh, eid, nbase_func, conn)

      do i = 1, nsurf
        ckey = get_key_surf(nbase_func, i, conn)
        is_exist = .false.
        call monolis_hash_get(hash_tree, ckey, in, is_exist)
        if(is_exist)then
          in = in + 1
        else
          in = 1
        endif
        call monolis_hash_push(hash_tree, ckey, in, is_pushed, is_exist)
      enddo
    enddo

    allocate(is_inner(nsurf,mesh%nelem), source = 0)
    if(present(is_open_surf)) allocate(is_open_surf(nsurf,mesh%nelem), source = 1)

    do eid = 1, mesh%nelem
      call monolis_get_connectivity(mesh, eid, nbase_func, conn)

      do i = 1, nsurf
        ckey = get_key_surf(nbase_func, i, conn)
        is_exist = .false.
        call monolis_hash_get(hash_tree, ckey, in, is_exist)
        if(.not. is_exist) stop "error: monolis_get_surf_node_tet"
        if(in == 2)then
          is_inner(i,eid) = 1
          if(present(is_open_surf)) is_open_surf(i,eid) = 0
        endif
      enddo
    enddo

    call monolis_hash_finalize(hash_tree)

    allocate(is_surf_node(mesh%nnode), source = 0)

    do eid = 1, mesh%nelem
      call monolis_get_connectivity(mesh, eid, nbase_func, conn)
      do i = 1, nsurf
        if(is_inner(i,eid) == 0)then
          do j = 1, nsurf_node
            if(nbase_func == 4)then
              in = conn(i341(j,i))
            elseif(nbase_func == 8)then
              in = conn(i361(j,i))
            endif
            is_surf_node(in) = 1
          enddo
        endif
      enddo
    enddo
  end subroutine monolis_get_surf_node

  function get_key_surf(nbase_func, i, conn)
    implicit none
    integer(kint) :: nbase_func, i, conn(:), array(4)
    integer(kint) :: i1, i2, i3, i4
    character :: c1*9, c2*9, c3*9, get_key_surf*27

    if(nbase_func == 4)then
      i1 = conn(i341(1,i))
      i2 = conn(i341(2,i))
      i3 = conn(i341(3,i))
      i4 = 2100000000
    elseif(nbase_func == 8)then
      i1 = conn(i361(1,i))
      i2 = conn(i361(2,i))
      i3 = conn(i361(3,i))
      i4 = conn(i361(4,i))
    else
      stop "error nbase_func"
    endif

    array(1) = i1
    array(2) = i2
    array(3) = i3
    array(4) = i4
    call monolis_qsort_int(array, 1, 4)

    write(c1,"(i9.9)")array(1)
    write(c2,"(i9.9)")array(2)
    write(c3,"(i9.9)")array(3)
    get_key_surf = c1//c2//c3
  end function get_key_surf
end module mod_monolis_dbc_all_util
