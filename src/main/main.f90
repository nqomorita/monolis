program main
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_solve
  use mod_monolis_hash
  implicit none
  type(monolis_prm) :: monoPRM
  type(monolis_com) :: monoCOM
  type(monolis_mat) :: monoMAT
  integer(kind=kint) :: key, val
  logical :: is_pushed, is_exist

  !call monolis_initialize(monoPRM, monoCOM, monoMAT)
  call monolis_hash_init()

  key = 321
  val = 123
  call monolis_hash_push(key, val, is_pushed, is_exist)
  write(*,*)"is_pushed, is_exist", is_pushed, is_exist

  key = 321
  call monolis_hash_get(key, val, is_exist)
  write(*,*)"val, is_exist", val, is_exist

  key = 111
  call monolis_hash_get(key, val, is_exist)
  write(*,*)"val, is_exist", val, is_exist

  call monolis_hash_finalize()
  !call monolis_finalize(monoPRM, monoCOM, monoMAT)

end program main
