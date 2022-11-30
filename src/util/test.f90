program monolis_test
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_util_debug
  use mod_monolis_neighbor_search
  implicit none
  type(type_monolis_neighbor_search) :: monolis_nbsearch
  integer(kint) :: div(3), nngrp, eid, i
  real(kdouble) :: BB(6), pos(3), ths
  integer(kint), allocatable :: ngrp(:)

  call monolis_global_initialize()

  call monolis_set_debug(.true.)
  call monolis_debug_header("test section")

  BB(1) = 0.0d0; BB(2) = 3.0d0
  BB(3) = 0.0d0; BB(4) = 3.0d0
  BB(5) = 0.0d0; BB(6) = 3.0d0
  div(1) = 3
  div(2) = 3
  div(3) = 3
  call monolis_neighbor_search_init(monolis_nbsearch, BB, div)

  call monolis_debug_header("monolis_neighbor_search_push")
  ths = 0.1d0

  eid = 101
  BB(1) = 0.0d0 + ths; BB(2) = 2.0d0 - ths
  BB(3) = 0.0d0 + ths; BB(4) = 2.0d0 - ths
  BB(5) = 0.0d0 + ths; BB(6) = 2.0d0 - ths
  call monolis_neighbor_search_push(monolis_nbsearch, BB, eid)

  eid = 102
  BB(1) = 1.0d0 + ths; BB(2) = 3.0d0 - ths
  BB(3) = 0.0d0 + ths; BB(4) = 2.0d0 - ths
  BB(5) = 0.0d0 + ths; BB(6) = 2.0d0 - ths
  call monolis_neighbor_search_push(monolis_nbsearch, BB, eid)

  eid = 103
  BB(1) = 1.0d0 + ths; BB(2) = 3.0d0 - ths
  BB(3) = 1.0d0 + ths; BB(4) = 3.0d0 - ths
  BB(5) = 0.0d0 + ths; BB(6) = 2.0d0 - ths
  call monolis_neighbor_search_push(monolis_nbsearch, BB, eid)

  eid = 104
  BB(1) = 0.0d0 + ths; BB(2) = 2.0d0 - ths
  BB(3) = 1.0d0 + ths; BB(4) = 3.0d0 - ths
  BB(5) = 0.0d0 + ths; BB(6) = 2.0d0 - ths
  call monolis_neighbor_search_push(monolis_nbsearch, BB, eid)

  eid = 105
  BB(1) = 0.0d0 + ths; BB(2) = 2.0d0 - ths
  BB(3) = 0.0d0 + ths; BB(4) = 2.0d0 - ths
  BB(5) = 1.0d0 + ths; BB(6) = 3.0d0 - ths
  call monolis_neighbor_search_push(monolis_nbsearch, BB, eid)

  eid = 106
  BB(1) = 1.0d0 + ths; BB(2) = 3.0d0 - ths
  BB(3) = 0.0d0 + ths; BB(4) = 2.0d0 - ths
  BB(5) = 1.0d0 + ths; BB(6) = 3.0d0 - ths
  call monolis_neighbor_search_push(monolis_nbsearch, BB, eid)

  eid = 107
  BB(1) = 1.0d0 + ths; BB(2) = 3.0d0 - ths
  BB(3) = 1.0d0 + ths; BB(4) = 3.0d0 - ths
  BB(5) = 1.0d0 + ths; BB(6) = 3.0d0 - ths
  call monolis_neighbor_search_push(monolis_nbsearch, BB, eid)

  eid = 108
  BB(1) = 0.0d0 + ths; BB(2) = 2.0d0 - ths
  BB(3) = 1.0d0 + ths; BB(4) = 3.0d0 - ths
  BB(5) = 1.0d0 + ths; BB(6) = 3.0d0 - ths
  call monolis_neighbor_search_push(monolis_nbsearch, BB, eid)

  !do i = 1, div(1)*div(2)*div(3)
  !  write(*,*) monolis_nbsearch%cell(i)%nid
  !  if(monolis_nbsearch%cell(i)%nid /= 0) write(*,*) monolis_nbsearch%cell(i)%id
  !  write(*,*)monolis_nbsearch%cell(i)%nid
  !  write(*,*)monolis_nbsearch%cell(i)%id
  !enddo

  call monolis_debug_header("monolis_neighbor_search_get_by_position")
  pos(1) = 0.0d0; pos(2) = 0.0d0; pos(3) = 0.0d0
  call monolis_neighbor_search_get_by_position(monolis_nbsearch, pos, nngrp, ngrp)
  write(*,*)nngrp, ngrp
  deallocate(ngrp)

  BB(1) = 0.0d0; BB(2) = 1.0d0
  BB(3) = 0.0d0; BB(4) = 1.0d0
  BB(5) = 0.0d0; BB(6) = 1.0d0
  call monolis_neighbor_search_get_by_bb(monolis_nbsearch, BB, nngrp, ngrp)
  write(*,*)nngrp, ngrp

  call monolis_neighbor_search_finalize(monolis_nbsearch)

  call monolis_global_finalize()

end program monolis_test
