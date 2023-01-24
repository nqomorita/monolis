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
  logical, save :: is_factored = .false.
  logical, save :: is_self = .false.

#ifdef WITH_MUMPS
  type (dmumps_struc), save :: mumps
#endif

contains

  subroutine monolis_precond_MUMPS_setup_local()
   implicit none
   is_self = .true.
  end subroutine monolis_precond_MUMPS_setup_local

  subroutine monolis_precond_mumps_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: NDOF, NZ

#ifdef WITH_MUMPS
    if(is_factored) return
    is_factored = .true.
    !is_self = .true.

    !> initialize
    mumps%JOB = -1
    mumps%COMM = monoCOM%comm
    if(is_self) mumps%COMM = mpi_comm_self
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
    mumps%N = NDOF*monoMAT%N

    call monolis_precond_mumps_get_offset(mumps%N)
    call monolis_precond_mumps_get_nz(monoMAT, NZ, is_self)

    allocate(mumps%IRN_loc(NDOF*NDOF*NZ))
    allocate(mumps%JCN_loc(NDOF*NDOF*NZ))
    allocate(mumps%A(NDOF*NDOF*NZ), source = 0.0d0)
    call monolis_precond_mumps_get_loc(monoMAT, monoCOM, &
      & mumps%IRN_loc, mumps%JCN_loc, mumps%A, is_self)

    mumps%JOB = 4
    mumps%NZ = NDOF*NDOF*NZ
    mumps%NZ_loc = NDOF*NDOF*NZ
    mumps%A_loc => mumps%A

    if(.not. is_self)then
      call monolis_allreduce_I1(mumps%N,  monolis_sum, monolis_global_comm())
      call monolis_allreduce_I1(mumps%NZ, monolis_sum, monolis_global_comm())
    endif

    if(monolis_global_myrank() == 0 .or. is_self)then
      allocate(mumps%RHS(mumps%N), source = 0.0d0)
    else
      allocate(mumps%RHS(1), source = 0.0d0)
    endif

    !> Output log level
    mumps%ICNTL(4) = 0
    mumps%ICNTL(14) = 80
    mumps%ICNTL(18) = 3

    call DMUMPS(mumps)
#else
    !call monolis_warning_header("monolis_precond_mumps_setup: MUMPS is NOT enabled")
    write(*,*)"* monolis_precond_mumps_setup: MUMPS is NOT enabled"
    stop
#endif
  end subroutine monolis_precond_mumps_setup

  subroutine monolis_precond_mumps_apply(monoPRM, monoCOM, monoMAT, X, Y)
    use mod_monolis_linalg_util
    use mod_monolis_matvec
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i
    real(kdouble) :: X(:), Y(:), t1, t2

#ifdef WITH_MUMPS
    !> Output log level
    mumps%ICNTL(4) = 0
    !> solution
    mumps%JOB = 3
    mumps%ICNTL(18) = 3

!Y = 1.0d0
!call monolis_matvec(monoCOM, monoMAT, Y, X, t1, t2)
!Y = 0.0d0

    call monolis_precond_mumps_set_rhs(monoMAT%NDOF*monoMAT%N, X, mumps%RHS)
    call DMUMPS(mumps)
    call monolis_precond_mumps_get_rhs(monoMAT%NDOF*monoMAT%N, mumps%RHS, Y)
    call monolis_update_R(monoCOM, monoMAT%NDOF, Y, t1)

!write(*,"(1p4e12.5)")Y
!stop
#endif
  end subroutine monolis_precond_mumps_apply

  subroutine monolis_precond_mumps_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%is_prec_stored) return

    is_factored = .false.

#ifdef WITH_MUMPS
    if(allocated(offset_list)) deallocate(offset_list)
    if(allocated(offset_counts)) deallocate(offset_counts)

    if(associated(mumps%IRN_loc)) deallocate(mumps%IRN_loc)
    if(associated(mumps%JCN_loc)) deallocate(mumps%JCN_loc)
    if(associated(mumps%A)) deallocate(mumps%A)
    if(associated(mumps%RHS)) deallocate(mumps%RHS)
#endif
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
    use mod_monolis_linalg_util
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_com) :: monoCOM
    integer(kint) :: i, j, in, jS, jE
    integer(kint) :: idof, jdof, jn, kn, NDOF, nprocs, myrank, i2, jn2, offset
    integer(kint), allocatable :: offset_id(:), offset_rank(:)
    integer(kint), pointer :: IRN(:), JCN(:)
    real(kdouble) :: A(:)
    logical :: is_self

    NDOF = monoMAT%ndof
    nprocs = monolis_global_commsize()
    myrank = monolis_global_myrank()

    offset = offset_list(monolis_global_myrank()+1)
    allocate(offset_id(monoMAT%NP), source = 0)
    allocate(offset_rank(monoMAT%NP), source = -1)

    do i = 1, monoMAT%NP
      offset_id(i) = i
    enddo
    call monolis_update_I(monoCOM, 1, offset_id)

    do i = 1, monoCOM%recv_n_neib
      jS = monoCOM%recv_index(i-1)+1
      jE = monoCOM%recv_index(i)
      do j = jS, jE
        in = monoCOM%recv_item(j)
        offset_rank(in) = monoCOM%recv_neib_pe(i)
      enddo
    enddo

    do i = 1, monoMAT%NP
      if(offset_rank(i) == -1)then
        offset_id(i) = offset_id(i) + offset_list(myrank+1)/monoMAT%NDOF
      else
        offset_id(i) = offset_id(i) + offset_list(offset_rank(i)+1)/monoMAT%NDOF
      endif
    enddo

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
            i2 = offset_id(i)
            jn2 = offset_id(jn)
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
    if(monolis_global_commsize() == 1 .or. is_self)then
      RHS = X(1:N)
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
    if(monolis_global_commsize() == 1 .or. is_self)then
      Y(1:N) = RHS
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