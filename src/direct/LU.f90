module mod_monolis_direct
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  implicit none

  integer(kind=kint), pointer :: perm(:)  => null()
  integer(kind=kint), pointer :: iperm(:) => null()
  logical, save :: isEntire = .false.

contains

  subroutine monolis_solver_LU(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMATF
    real(kind=kdouble) :: t1, t2, tset, tsol, tcomm, R2

    !allocate( perm(NP))
    !allocate(iperm(NP))

    if(monoCOM%myrank == 0) write(*,"(a)")"   ** analysis phase"
    if(monoCOM%commsize == 1) isEntire = .true.

    t1 = monolis_wtime()
    !call monolis_get_fillin(monoMAT, hecT, indexU, itemU, NPU)
    t2 = monolis_wtime()
    tset = tset + t2 - t1
    if(monoCOM%myrank == 0) write(*,"(a,1pe11.4)")"    * fill-in    time: ", t2-t1

    t1 = monolis_wtime()
    !call monolis_matrix_copy_with_fillin(hecMESH, monoMATF, indexU, itemU, AU, NPU)
    t2 = monolis_wtime()
    tset = tset + t2 - t1
    if(monoCOM%myrank == 0) write(*,"(a,1pe11.4)")"    * reallocate time: ", t2-t1

  end subroutine monolis_solver_LU
end module mod_monolis_direct
