module monolis
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_hecmw
end module monolis

subroutine monolis_solve_hecmw()
  use monolis
  implicit none
  type(monolis_prm) :: monoPRM
  type(monolis_com) :: monoCOM
  type(monolis_mat) :: monoMAT

  call monolis_solve_hecmw_inner()

end subroutine monolis_solve_hecmw
