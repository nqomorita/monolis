module mod_monolis_neighbor_search
  use mod_monolis_prm
  use mod_monolis_stdlib

  private

  public :: monolis_neighbor_search_init
  public :: monolis_neighbor_search_finalize
  public :: monolis_neighbor_search_get
  public :: monolis_neighbor_search_push
  public :: type_monolis_neighbor_search

  type type_monolis_neighbor_search_main
    integer(kint) :: nid = 0
    integer(kint), pointer :: id(:) => null()
  end type type_monolis_neighbor_search_main

  type type_monolis_neighbor_search
    integer(kint) :: div(3)
    real(kdouble) :: BB(6) !> xmin, xmax, ymin, ymax, zmin, zmax
    real(kdouble) :: dx(3) !> dx, dy, dz
    type(type_monolis_neighbor_search_main), allocatable :: cell(:)
  end type type_monolis_neighbor_search

  real(kdouble), parameter :: ths = 1.0d-3

contains

  subroutine monolis_neighbor_search_init(monolis_nbsearch, BB, div)
    implicit none
    type(type_monolis_neighbor_search) :: monolis_nbsearch
    integer(kint) :: div(3)
    real(kdouble) :: BB(6)
    monolis_nbsearch%div = div
    monolis_nbsearch%BB = BB
    monolis_nbsearch%dx(1) = (BB(2)-BB(1))/dble(div(1))
    monolis_nbsearch%dx(2) = (BB(4)-BB(3))/dble(div(2))
    monolis_nbsearch%dx(3) = (BB(6)-BB(5))/dble(div(3))
    allocate(monolis_nbsearch%cell(div(1)*div(2)*div(3)))
  end subroutine monolis_neighbor_search_init

  subroutine monolis_neighbor_search_push(monolis_nbsearch, BB, id)
    implicit none
    type(type_monolis_neighbor_search) :: monolis_nbsearch
    integer(kint) :: id, x, y, z, in, imin(3), imax(3)
    real(kdouble) :: BB(6), pos(3)

    pos(1) = BB(1) - ths
    pos(2) = BB(3) - ths
    pos(3) = BB(5) - ths
    call get_int_coordinate(monolis_nbsearch, pos, imin)

    pos(1) = BB(2) + ths
    pos(2) = BB(4) + ths
    pos(3) = BB(6) + ths
    call get_int_coordinate(monolis_nbsearch, pos, imax)

    do z = 1, imax(1) - imin(1) + 1
    do y = 1, imax(2) - imin(2) + 1
    do x = 1, imax(3) - imin(3) + 1
      in = 1
      call monolis_neighbor_search_push_main(monolis_nbsearch, in, id)
    enddo
    enddo
    enddo
  end subroutine monolis_neighbor_search_push

  subroutine monolis_neighbor_search_push_main(monolis_nbsearch, in, id)
    implicit none
    type(type_monolis_neighbor_search) :: monolis_nbsearch
    integer(kint) :: in, id, add(1), nid
    nid = monolis_nbsearch%cell(in)%nid
    add = id
    call monolis_pointer_reallocate_integer(monolis_nbsearch%cell(in)%id, nid, 1, add)
  end subroutine monolis_neighbor_search_push_main

  subroutine monolis_neighbor_search_get(monolis_nbsearch, pos, nid, id)
    implicit none
    type(type_monolis_neighbor_search) :: monolis_nbsearch
    integer(kint) :: nid, morton_id
    integer(kint), pointer :: id(:)
    real(kdouble) :: pos(3)
    logical :: is_in

    call BB_check(monolis_nbsearch, pos, is_in)
    if(.not. is_in)then
      nid = 0
      nullify(id)
      return
    endif

    call get_z_index_by_position(monolis_nbsearch, pos, morton_id)

    nid = monolis_nbsearch%cell(morton_id)%nid
    if(nid < 1)then
      nullify(id)
    else
      id => monolis_nbsearch%cell(morton_id)%id
    endif
  end subroutine monolis_neighbor_search_get

  subroutine monolis_neighbor_search_finalize(monolis_nbsearch)
    implicit none
    type(type_monolis_neighbor_search) :: monolis_nbsearch
    deallocate(monolis_nbsearch%cell)
  end subroutine monolis_neighbor_search_finalize

  subroutine BB_check(monolis_nbsearch, pos, is_in)
    implicit none
    type(type_monolis_neighbor_search) :: monolis_nbsearch
    real(kdouble) :: pos(3), BB(6)
    logical :: is_in

    is_in = .false.
    BB = monolis_nbsearch%BB
    if( BB(1)-ths <= pos(1) .and. pos(1) <= BB(2)+ths .and. &
      & BB(3)-ths <= pos(2) .and. pos(2) <= BB(4)+ths .and. &
      & BB(5)-ths <= pos(3) .and. pos(3) <= BB(6)+ths )then
      is_in = .true.
    endif
  end subroutine BB_check

  subroutine get_z_index_by_position(monolis_nbsearch, pos, morton_id)
    implicit none
    type(type_monolis_neighbor_search) :: monolis_nbsearch
    integer(kint) :: morton_id, id(3), div(3)
    real(kdouble) :: pos(3)

    call get_int_coordinate(monolis_nbsearch, pos, id)
    div = monolis_nbsearch%div
    morton_id = id(3)*div(1)*div(2) + id(2)*div(1) + id(1)
  end subroutine get_z_index_by_position

  subroutine get_int_coordinate(monolis_nbsearch, pos, id)
    implicit none
    type(type_monolis_neighbor_search) :: monolis_nbsearch
    integer(kint) :: id(3)
    real(kdouble) :: pos(3)

    id(1) = int((pos(1) - monolis_nbsearch%BB(1))/monolis_nbsearch%dx(1))
    id(2) = int((pos(2) - monolis_nbsearch%BB(3))/monolis_nbsearch%dx(2))
    id(3) = int((pos(3) - monolis_nbsearch%BB(5))/monolis_nbsearch%dx(3))
    if(id(1) < 0) id(1) = 0
    if(id(2) < 0) id(2) = 0
    if(id(3) < 0) id(3) = 0
  end subroutine get_int_coordinate
end module mod_monolis_neighbor_search
