module mod_monolis_p_refiner
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  implicit none

contains

  subroutine monolis_p_refine(mesh, mesh_ref)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh, mesh_ref

  end subroutine monolis_p_refine

end module mod_monolis_p_refiner
