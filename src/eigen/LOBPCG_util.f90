module blopex_fortran_hold_vars
  use mod_monolis_utils
  use mod_monolis_def_solver
  use mod_monolis_def_mat
  use iso_c_binding
  integer(c_int) :: N_hold, M_hold
  type(monolis_prm), pointer :: monoPRM_hold
  type(monolis_com), pointer :: monoCOM_hold
  type(monolis_mat), pointer :: monoMAT_hold
  real(kdouble), pointer :: mass_hold(:)
  real(kdouble), pointer :: filter_hold(:)
end module blopex_fortran_hold_vars

subroutine blopex_fortran_opA(dum, a, b)
  use blopex_fortran_hold_vars
  use mod_monolis_matvec
  use iso_c_binding
  implicit none
  integer(c_int) :: dum, i, j, jS, jE
  real(c_double) :: a(N_hold*M_hold), b(N_hold*M_hold), t(N_hold*M_hold)
  real(c_double) :: t1, t2

  do j = 1, M_hold
    jS = N_hold*(j-1) + 1
    jE = N_hold*(j)
    do i = 1, N_hold
      t(i) = a(jS-1+i)*filter_hold(i)
    enddo
    call monolis_matvec(monoCOM_hold, monoMAT_hold, t, b(jS:jE), t1, t2)
  enddo
end subroutine blopex_fortran_opA

subroutine blopex_fortran_opB(dum, a, b)
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none
  integer(c_int) :: dum, i, j, jS, jE
  real(c_double) :: a(N_hold*M_hold), b(N_hold*M_hold)

  do j = 1, M_hold
    jS = N_hold*(j-1) + 1
    jE = N_hold*(j)
    do i = jS, jE
      b(i) = mass_hold(i-jS+1)*a(i)*filter_hold(i-jS+1)
    enddo
  enddo
end subroutine blopex_fortran_opB

subroutine blopex_fortran_opT(dum, a, b)
  use blopex_fortran_hold_vars
  use iso_c_binding
  use mod_monolis_precond
  implicit none
  integer(c_int) :: dum, i, j, jS, jE
  real(c_double) :: a(N_hold*M_hold), b(N_hold*M_hold), w(N_hold), t(N_hold*M_hold), tcomm

  do j = 1, M_hold
    jS = N_hold*(j-1) + 1
    jE = N_hold*(j)
    do i = 1, N_hold
      t(i) = a(jS-1+i)*filter_hold(i)
    enddo
    call monolis_precond_apply(monoPRM_hold, monoCOM_hold, monoMAT_hold, t, b(jS:jE))
  enddo
end subroutine blopex_fortran_opT
