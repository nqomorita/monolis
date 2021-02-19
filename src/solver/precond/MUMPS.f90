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

  logical, save :: is_factored = .false.

#ifdef WITH_MUMPS
  type (dmumps_struc), save :: mumps
#endif

contains

  subroutine monolis_precond_mumps_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: NDOF, NZ

#ifdef WITH_MUMPS
    if(is_factored) return
    is_factored = .true.

    !> initialize
    mumps%JOB = -1
    mumps%COMM = mpi_comm_self !monoCOM%comm
    !mumps%SYM = mumps_mat_spd
    mumps%SYM = mumps_mat_asym
    !> parallel fatorization, 0:serial, 1:parallel
    mumps%PAR = 1
    !> ordering: 0:auto, 1:seq, 2:par
    mumps%ICNTL(28) = 0
    !> seq ord: 0:AMD, 1:USER, 2:AMF, 3:scotch, 4:pord, 5:metis, 6:QAMD, 7:auto
    mumps%ICNTL(7) = 7
    !> par ord: 0:auto, 1:ptscotch, 2:parmetis
    mumps%ICNTL(29) = 0
    !> num of OMP
    mumps%ICNTL(16) = 1
    !> iterative refinement
    mumps%ICNTL(10) = 2
    mumps%CNTL(2) = 1.0e-8
    !> Out-Of-Core: 0:IN-CORE only, 1:OOC
    mumps%ICNTL(22) = 0
    !> Distributed assembled matrix input
    mumps%ICNTL(18) = 3
    mumps%ICNTL(14) = 60

    call DMUMPS(mumps)

    !> factorization
    NDOF = monoMAT%ndof
    call monolis_precond_mumps_get_nz(monoMAT, NZ)

    allocate(mumps%IRN_loc(NDOF*NDOF*NZ))
    allocate(mumps%JCN_loc(NDOF*NDOF*NZ))
    allocate(mumps%IRN(NDOF*NDOF*NZ))
    allocate(mumps%JCN(NDOF*NDOF*NZ))
    allocate(mumps%A(NDOF*NDOF*NZ), source = 0.0d0)
    allocate(mumps%RHS(NDOF*monoMAT%N), source = 0.0d0)
    call monolis_precond_mumps_get_loc(monoMAT, mumps%IRN_loc, mumps%JCN_loc, mumps%A)

    mumps%IRN = mumps%IRN_loc
    mumps%JCN = mumps%JCN_loc

    mumps%JOB = 4
    mumps%N = NDOF*monoMAT%N
    mumps%NZ = NDOF*NDOF*NZ
    mumps%NZ_loc = NDOF*NDOF*NZ

    !> Output log level
    mumps%ICNTL(4) = 0
    mumps%ICNTL(14) = 80

    call DMUMPS(mumps)
#endif
  end subroutine monolis_precond_mumps_setup

  subroutine monolis_precond_mumps_apply(monoPRM, monoCOM, monoMAT, X, Y)
    use mod_monolis_linalg_util
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i
    real(kdouble) :: X(:), Y(:), tcomm

#ifdef WITH_MUMPS
    !> Output log level
    mumps%ICNTL(4) = 0
    !> solution
    mumps%JOB = 3

    do i = 1, mumps%N
      mumps%RHS(i) = X(i)
    enddo

    call DMUMPS(mumps)

    do i = 1, mumps%N
      Y(i) = mumps%RHS(i)
    enddo
#endif
  end subroutine monolis_precond_mumps_apply

  subroutine monolis_precond_mumps_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

  end subroutine monolis_precond_mumps_clear

  subroutine monolis_precond_mumps_get_nz(monoMAT, NZ)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, in, jS, jE, NZ
    integer(kint) :: idof, jdof, jn, kn, NDOF

    NDOF = monoMAT%ndof
    NZ = 0

    in = 0
    aa:do i = 1, monoMAT%N
      jS = monoMAT%index(i-1) + 1
      jE = monoMAT%index(i)
      do j = jS, jE
        jn = monoMAT%item(j)
        if(monoMAT%N < jn) cycle aa
        NZ = NZ + 1
      enddo
    enddo aa
  end subroutine monolis_precond_mumps_get_nz

  subroutine monolis_precond_mumps_get_loc(monoMAT, IRN, JCN, A)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, in, jS, jE
    integer(kint) :: idof, jdof, jn, kn, NDOF
    integer(kint), pointer :: IRN(:), JCN(:)
    real(kdouble) :: A(:)

    NDOF = monoMAT%ndof

    in = 0
    aa:do i = 1, monoMAT%N
      jS = monoMAT%index(i-1) + 1
      jE = monoMAT%index(i)
      do j = jS, jE
        jn = monoMAT%item(j)
        if(monoMAT%N < jn) cycle aa
        do idof = 1, NDOF
        do jdof = 1, NDOF
          kn = NDOF*NDOF*(j-1) + NDOF*(idof-1) + jdof
          in = in + 1
          A(in) = monoMAT%A(kn)
          IRN(in) = NDOF*(i  - 1) + idof
          JCN(in) = NDOF*(jn - 1) + jdof
        enddo
        enddo
      enddo
    enddo aa
  end subroutine monolis_precond_mumps_get_loc
end module mod_monolis_precond_mumps