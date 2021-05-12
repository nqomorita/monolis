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

  integer(kint) :: i341(3,4), i361(4,6)
  type(type_monolis_hash_tree) :: hash_tree

contains

  subroutine monolis_get_surf_node(mesh, nbase_func, nsurf, nsurf_node, is_surf_node)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh
    integer(kint) :: nbase_func, conn(nbase_func), nsurf, nsurf_node
    integer(kint) :: i, j, in, eid, tmp
    character :: ckey*27
    logical :: is_exist, is_pushed
    integer(kint), allocatable :: is_inner(:,:), is_surf_node(:)

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

    do eid = 1, mesh%nelem
      call monolis_get_connectivity(mesh, eid, nbase_func, conn)

      do i = 1, nsurf
        ckey = get_key_surf(nbase_func, i, conn)
        is_exist = .false.
        call monolis_hash_get(hash_tree, ckey, in, is_exist)
        if(.not. is_exist) stop "error: monolis_get_surf_node_tet"
        if(in == 2) is_inner(i,eid) = 1
      enddo
    enddo

    call monolis_hash_finalize(hash_tree)

    allocate(is_surf_node(mesh%nnode), source = 0)

    do eid = 1, mesh%nelem
      call monolis_get_connectivity(mesh, eid, nbase_func, conn)
      do i = 1, nsurf
        if(is_inner(i,eid) == 1)then
          do j = 1, nsurf_node
            in = conn(i341(j,i))
            is_surf_node(in) = 1
          enddo
        endif
      enddo
    enddo
  end subroutine monolis_get_surf_node

  function get_key_surf(nbase_func, i, conn)
    implicit none
    integer(kint) :: nbase_func, i, conn(:), array(3)
    integer(kint) :: i1, i2, i3, i4
    character :: c1*9, c2*9, c3*9, get_key_surf*27

    if(nbase_func == 4)then
      i1 = conn(i341(1,i))
      i2 = conn(i341(2,i))
      i3 = conn(i341(3,i))
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
    call monolis_qsort_int(array, 1, 3)

    write(c1,"(i9.9)")array(1)
    write(c2,"(i9.9)")array(2)
    write(c3,"(i9.9)")array(3)
    get_key_surf = c1//c2//c3
  end function get_key_surf
end module mod_monolis_dbc_all_util
