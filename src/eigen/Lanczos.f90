module mod_monolis_eigen_lanczos
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

contains

  subroutine monolis_eigen_lanczos(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_lanczos")

  end subroutine monolis_eigen_lanczos

end module mod_monolis_eigen_lanczos
