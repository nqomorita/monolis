module mod_monolis_def_solver_util
  use mod_monolis_utils
  use mod_monolis_def_struc

  implicit none

  private
  public :: monolis_param_set_method
  public :: monolis_param_set_precond
  public :: monolis_param_set_maxiter
  public :: monolis_param_set_tolerance
  public :: monolis_param_set_scaling
  public :: monolis_param_set_reordering
  public :: monolis_param_set_init_x
  public :: monolis_param_set_sym_matrix
  public :: monolis_param_set_debug
  public :: monolis_param_set_check_diag
  public :: monolis_param_set_show_iterlog
  public :: monolis_param_set_show_time
  public :: monolis_param_set_show_summary

contains

  subroutine monolis_param_set_method(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: param
    monolis%PRM%method = param
  end subroutine monolis_param_set_method

  subroutine monolis_param_set_precond(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: param
    monolis%PRM%precond = param
  end subroutine monolis_param_set_precond

  subroutine monolis_param_set_maxiter(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: param
    monolis%PRM%maxiter = param
  end subroutine monolis_param_set_maxiter

  subroutine monolis_param_set_tolerance(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    real(kdouble) :: param
    monolis%PRM%tol = param
  end subroutine monolis_param_set_tolerance

  subroutine monolis_param_set_scaling(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_scaling = param
  end subroutine monolis_param_set_scaling

  subroutine monolis_param_set_reordering(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_reordering = param
  end subroutine monolis_param_set_reordering

  subroutine monolis_param_set_init_x(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_init_x = param
  end subroutine monolis_param_set_init_x

  subroutine monolis_param_set_sym_matrix(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_sym_matrix = param
  end subroutine monolis_param_set_sym_matrix

  subroutine monolis_param_set_debug(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_debug = param
  end subroutine monolis_param_set_debug

  subroutine monolis_param_set_check_diag(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_check_diag = param
  end subroutine monolis_param_set_check_diag

  subroutine monolis_set_performance_measurement(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_measurement = param
  end subroutine monolis_set_performance_measurement

  subroutine monolis_param_set_show_iterlog(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%show_iterlog = param
  end subroutine monolis_param_set_show_iterlog

  subroutine monolis_param_set_show_time(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%show_time = param
  end subroutine monolis_param_set_show_time

  subroutine monolis_param_set_show_summary(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%show_summary = param
  end subroutine monolis_param_set_show_summary
end module mod_monolis_def_solver_util
