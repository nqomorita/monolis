module mod_monolis_precond_mumps
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg_com

  implicit none

#ifdef WITH_MUMPS
  include 'dmumps_struc.h'
#endif

  integer(kint), parameter :: mumps_mat_csr = 1
  integer(kint), parameter :: mumps_mat_coo = 2
  integer(kint), parameter :: mumps_mat_asym = 0
  integer(kint), parameter :: mumps_mat_spd = 1
  integer(kint), parameter :: mumps_mat_sym = 2

  integer(kint), save, allocatable :: offset_list(:)
  integer(kint), save, allocatable :: offset_counts(:)
  !integer(kint), save :: offset
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
    mumps%COMM = monoCOM%comm
    !mumps%COMM = mpi_comm_self
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
    call monolis_precond_mumps_get_nz(monoMAT, NZ, is_self = .false.)

    allocate(mumps%IRN_loc(NDOF*NDOF*NZ))
    allocate(mumps%JCN_loc(NDOF*NDOF*NZ))
    allocate(mumps%A(NDOF*NDOF*NZ), source = 0.0d0)
    call monolis_precond_mumps_get_loc(monoMAT, monoCOM, &
      & mumps%IRN_loc, mumps%JCN_loc, mumps%A, is_self = .false.)

    mumps%JOB = 4
    mumps%N = NDOF*monoMAT%N
    mumps%NZ = NDOF*NDOF*NZ
    mumps%NZ_loc = NDOF*NDOF*NZ
    mumps%A_loc => mumps%A

    call monolis_precond_mumps_get_offset(mumps%N)
    call monolis_allreduce_I1(mumps%N,  monolis_sum, monolis_global_comm())
    call monolis_allreduce_I1(mumps%NZ, monolis_sum, monolis_global_comm())

    if(monolis_global_myrank() == 0)then
      allocate(mumps%RHS(mumps%N), source = 0.0d0)
    else
      allocate(mumps%RHS(1), source = 0.0d0)
    endif

    !> Output log level
    mumps%ICNTL(4) = 0
    mumps%ICNTL(14) = 80
    mumps%ICNTL(18) = 3

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
    mumps%ICNTL(18) = 3

    call monolis_precond_mumps_set_rhs(monoMAT%NDOF*monoMAT%N, X, mumps%RHS)
    call DMUMPS(mumps)
    call monolis_precond_mumps_get_rhs(monoMAT%NDOF*monoMAT%N, mumps%RHS, Y)
#endif
  end subroutine monolis_precond_mumps_apply

  subroutine monolis_precond_mumps_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
  end subroutine monolis_precond_mumps_clear

  subroutine monolis_precond_mumps_get_nz(monoMAT, NZ, is_self)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, in, jS, jE, NZ
    integer(kint) :: idof, jdof, jn, kn, NDOF
    logical :: is_self

    NDOF = monoMAT%ndof
    NZ = 0

    in = 0
    aa:do i = 1, monoMAT%N
      jS = monoMAT%index(i-1) + 1
      jE = monoMAT%index(i)
      do j = jS, jE
        jn = monoMAT%item(j)
        if(is_self .and. monoMAT%N < jn) cycle aa
        NZ = NZ + 1
      enddo
    enddo aa
  end subroutine monolis_precond_mumps_get_nz

  subroutine monolis_precond_mumps_get_loc(monoMAT, monoCOM, IRN, JCN, A, is_self)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_com) :: monoCOM
    integer(kint) :: i, j, in, jS, jE
    integer(kint) :: idof, jdof, jn, kn, NDOF, nprocs, i2, jn2
    integer(kint), pointer :: IRN(:), JCN(:)
    real(kdouble) :: A(:)
    logical :: is_self

    NDOF = monoMAT%ndof
    nprocs = monolis_global_commsize()

    in = 0
    aa:do i = 1, monoMAT%N
      jS = monoMAT%index(i-1) + 1
      jE = monoMAT%index(i)
      do j = jS, jE
        jn = monoMAT%item(j)
        if(is_self .and. monoMAT%N < jn) cycle aa
        do idof = 1, NDOF
        do jdof = 1, NDOF
          kn = NDOF*NDOF*(j-1) + NDOF*(idof-1) + jdof
          in = in + 1
          A(in) = monoMAT%A(kn)
          if(is_self .or. nprocs == 1)then
            IRN(in) = NDOF*(i  - 1) + idof
            JCN(in) = NDOF*(jn - 1) + jdof
          else
            i2 = monoCOM%global_node_id(i)
            jn2 = monoCOM%global_node_id(jn)
            IRN(in) = NDOF*(i2  - 1) + idof
            JCN(in) = NDOF*(jn2 - 1) + jdof
          endif
        enddo
        enddo
      enddo
    enddo aa
  end subroutine monolis_precond_mumps_get_loc

  subroutine monolis_precond_mumps_set_rhs(N, X, RHS)
    implicit none
    integer(kint) :: N
    real(kdouble) :: X(:), RHS(:)
    if(monolis_global_commsize() == 1)then
      RHS = X
    else
      call monolis_gatherv_R(X, N, &
        RHS, offset_counts, offset_list, &
        0, monolis_global_comm())
    endif
  end subroutine monolis_precond_mumps_set_rhs

  subroutine monolis_precond_mumps_get_rhs(N, RHS, Y)
    implicit none
    integer(kint) :: N
    real(kdouble) :: Y(:), RHS(:)
    if(monolis_global_commsize() == 1)then
      Y = RHS
    else
      call monolis_scatterv_R( &
        RHS, offset_counts, offset_list, &
        Y, N, &
        0, monolis_global_comm())
    endif
  end subroutine monolis_precond_mumps_get_rhs

  subroutine monolis_precond_mumps_get_offset(N_loc)
    implicit none
    integer(kint) :: nprocs, myrank, i, N_loc

    nprocs = monolis_global_commsize()
    myrank = monolis_global_myrank()

    allocate(offset_list(nprocs), source = 0)
    allocate(offset_counts(nprocs), source = 0)

    if(1 < nprocs)then
      call monolis_allgather_I1(N_loc, offset_counts, monolis_global_comm())
    endif

    offset_list(1) = 0
    do i = 1, nprocs - 1
      offset_list(i+1) = offset_list(i) + offset_counts(i)
    enddo
  end subroutine monolis_precond_mumps_get_offset
end module mod_monolis_precond_mumps