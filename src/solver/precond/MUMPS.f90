module mod_monolis_precond_mumps
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

#ifdef WITH_MUMPS
  include 'dmumps_struc.h'
#endif

  integer(kint), parameter :: mumps_mat_csr = 1
  integer(kint), parameter :: mumps_mat_coo = 2
  integer(kint), parameter :: mumps_mat_asym = 0
  integer(kint), parameter :: mumps_mat_spd = 1
  integer(kint), parameter :: mumps_mat_sym = 2

#ifdef WITH_MUMPS
  type (dmumps_struc), save :: mumps
#endif

contains

  subroutine monolis_precond_mumps_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: NDOF
    integer(kint), pointer :: IRN(:), JCN(:), MAP(:)

#ifdef WITH_MUMPS
    !> initialize
    mumps%JOB = -1
    mumps%COMM = monoCOM%comm
    mumps%SYM = mumps_mat_asym
    !> parallel fatorization, 0:serial, 1:parallel
    mumps%PAR = 1
    !> ordering: 0:auto, 1:seq, 2:par
    mumps%ICNTL(28) = 0
    !> seq ord: 0:AMD, 1:USER, 2:AMF, 3:scotch, 4:pord, 5:metis, 6:QAMD, 7:auto
    mumps%ICNTL(7) = 7
    !> par ord: 0:auto, 1:ptscotch, 2:parmetis
    mumps%ICNTL(29) = 0
    !> relaxation parameter
    mumps%ICNTL(14) = 20
    !> iterative refinement
    mumps%ICNTL(10) = 3
    mumps%CNTL(2) = 1.0e-8
    !> Out-Of-Core: 0:IN-CORE only, 1:OOC
    mumps%ICNTL(22) = 0
    !> Distributed assembled matrix input
    !mumps%ICNTL(18) = 3
    call DMUMPS(mumps)

    !> factorization
    NDOF = monoMAT%ndof
    allocate(IRN(NDOF*NDOF*monoMAT%NZ))
    allocate(JCN(NDOF*NDOF*monoMAT%NZ))
    allocate(mumps%A(NDOF*NDOF*monoMAT%NZ))
    !allocate(MAP(monoMAT%NZ))
    allocate(mumps%RHS(monoMAT%N*NDOF))
    call monolis_precond_mumps_get_loc(monoMAT, IRN, JCN, mumps%A)

    mumps%JOB = 4
    mumps%N = monoMAT%N*monoMAT%NDOF
    mumps%NZ = monoMAT%NZ
   ! mumps%NNZ_loc = monoMAT%NZ
   ! mumps%IRN_loc => IRN
   ! mumps%JCN_loc => JCN
    mumps%IRN => IRN
    mumps%JCN => JCN
    !mumps%A_loc => monoMAT%A

    call DMUMPS(mumps)
#endif
  end subroutine monolis_precond_mumps_setup

  subroutine monolis_precond_mumps_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kdouble) :: X(:), Y(:)

#ifdef WITH_MUMPS
    !> solution
    mumps%JOB = 3
    mumps%RHS = X
    call DMUMPS(mumps)
    Y = mumps%RHS
#endif
  end subroutine monolis_precond_mumps_apply

  subroutine monolis_precond_mumps_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

  end subroutine monolis_precond_mumps_clear

  subroutine monolis_precond_mumps_get_loc(monoMAT, IRN, JCN, A)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, in, jS, jE
    integer(kint) :: idof, jdof, jn, kn, NDOF
    integer(kint), pointer :: IRN(:), JCN(:)
    real(kdouble) :: A(:)

    NDOF = monoMAT%NDOF

    in = 1
    do i = 1, monoMAT%NP
      do idof = 1, monoMAT%NDOF
        jS = monoMAT%index(i-1) + 1
        jE = monoMAT%index(i)
        do j = jS, jE
          jn = monoMAT%item(j)
          do jdof = 1, NDOF
            kn = NDOF*NDOF*(j-1) + NDOF*(idof-1) + jdof
            A(in) = monoMAT%A(kn)
            JCN(in) = NDOF*(jn - 1) + jdof
            IRN(in) = NDOF*(i  - 1) + idof
            in = in + 1
          enddo
        enddo
      enddo
    enddo
  end subroutine monolis_precond_mumps_get_loc
end module mod_monolis_precond_mumps