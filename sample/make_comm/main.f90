program main
  use mod_monolis
  use mod_monolis_util_com
  implicit none
  type(monolis_structure) :: mat !> 疎行列変数
  !integer(kint) :: nnode, nelem, nbase_func, i
  !real(kdouble) :: resid
  !integer(kint), allocatable :: elem(:,:)
  integer(kint), allocatable :: nid(:)
  !real(kdouble), allocatable :: node(:,:), coef(:), B(:), X(:)
  !character :: fname*128

  call monolis_global_initialize()

  allocate(nid(3), source = 0)
  if(monolis_global_myrank() == 0)then
    nid(1) = 1
    nid(2) = 2
    nid(3) = 3
  elseif(monolis_global_myrank() == 1)then
    nid(1) = 3
    nid(2) = 4
    nid(3) = 2
  endif

  call monolis_com_get_comm_table(mat, 2, 3, nid)

  call monolis_global_finalize()

end program main
