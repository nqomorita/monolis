module mod_monolis_precond_mumps
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

  integer(kint), parameter :: mumps_mat_csr = 1
  integer(kint), parameter :: mumps_mat_coo = 2

  integer(kint), parameter :: mumps_mat_asym = 0
  integer(kint), parameter :: mumps_mat_spd = 1
  integer(kint), parameter :: mumps_mat_sym = 2

contains

  subroutine monolis_precond_MUMPS_setup_local(monoPREC)
    implicit none
    type(monolis_mat) :: monoPREC
    monoPREC%DMUMPS%is_self = .true.
  end subroutine monolis_precond_MUMPS_setup_local

  subroutine monolis_precond_mumps_setup(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_mat), target :: monoPREC
    integer(kint) :: NDOF, NZ
#ifdef WITH_MUMPS
    type (dmumps_struc), pointer :: mumps
    real(kdouble) :: t1, t2, t3, t4, t5

    t1 = monolis_get_time_global_sync()

    mumps => monoPREC%DMUMPS%mumps

    if(monoPREC%DMUMPS%is_factored) return
    monoPREC%DMUMPS%is_factored = .true.

    !> initialize
    mumps%JOB = -1
    mumps%COMM = monoCOM%comm
    if(monoPREC%DMUMPS%is_self) mumps%COMM = mpi_comm_self
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
    mumps%ICNTL(10) = 0
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

    call monolis_precond_mumps_get_offset(mumps%N, monoCOM%comm, &
      & monoPREC%DMUMPS%offset_list, monoPREC%DMUMPS%offset_counts)
    call monolis_precond_mumps_get_nz(monoMAT, NZ, monoPREC%DMUMPS%is_self)

    allocate(mumps%IRN_loc(NDOF*NDOF*NZ))
    allocate(mumps%JCN_loc(NDOF*NDOF*NZ))
    allocate(mumps%A(NDOF*NDOF*NZ), source = 0.0d0)
    call monolis_precond_mumps_get_loc(monoMAT, monoPREC, monoCOM, &
      & mumps%IRN_loc, mumps%JCN_loc, mumps%A, monoPREC%DMUMPS%is_self)

    mumps%JOB = 4
    mumps%NZ = NDOF*NDOF*NZ
    mumps%NZ_loc = NDOF*NDOF*NZ
    mumps%A_loc => mumps%A

    if(.not. monoPREC%DMUMPS%is_self)then
      call monolis_allreduce_I1(mumps%N,  monolis_mpi_sum, monoCOM%comm)
      call monolis_allreduce_I1(mumps%NZ, monolis_mpi_sum, monoCOM%comm)
    endif

    if(monolis_mpi_get_local_my_rank(monoCOM%comm) == 0 .or. monoPREC%DMUMPS%is_self)then
      allocate(mumps%RHS(mumps%N), source = 0.0d0)
    else
      allocate(mumps%RHS(1), source = 0.0d0)
    endif

    !> Output log level
    mumps%ICNTL(4) = 0
    mumps%ICNTL(14) = 80
    mumps%ICNTL(18) = 3

    call DMUMPS(mumps)

    t2 = monolis_get_time_global_sync()
    write(*,"(a,1pe10.3)")"analysis", t2 - t1
#else
    write(*,*)"* monolis_precond_mumps_setup: MUMPS is NOT enabled"
    stop
#endif
  end subroutine monolis_precond_mumps_setup

  subroutine monolis_precond_mumps_apply(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_mat), target :: monoPREC
    integer(kint) :: i
    real(kdouble) :: X(:), Y(:), t1, t2
#ifdef WITH_MUMPS
    type (dmumps_struc), pointer :: mumps

    mumps => monoPREC%DMUMPS%mumps

    !> Output log level
    mumps%ICNTL(4) = 0
    !> solution
    mumps%JOB = 3
    mumps%ICNTL(18) = 3

    call monolis_precond_mumps_set_rhs(monoMAT%NDOF*monoMAT%N, X, mumps%RHS, monoCOM%comm, &
      & monoPREC%DMUMPS%offset_list, monoPREC%DMUMPS%offset_counts, monoPREC%DMUMPS%is_self)
    call DMUMPS(mumps)
    call monolis_precond_mumps_get_rhs(monoMAT%NDOF*monoMAT%N, mumps%RHS, Y, monoCOM%comm, &
      & monoPREC%DMUMPS%offset_list, monoPREC%DMUMPS%offset_counts, monoPREC%DMUMPS%is_self)
    call monolis_mpi_update_R(monoCOM, monoMAT%NDOF, Y, t1)
#endif
  end subroutine monolis_precond_mumps_apply

  subroutine monolis_precond_mumps_clear(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_mat), target :: monoPREC
#ifdef WITH_MUMPS
    type(dmumps_struc), pointer :: mumps

    return
    !if(monoPRM%is_prec_stored) return

    mumps => monoPREC%DMUMPS%mumps

    monoPREC%DMUMPS%is_factored = .false.

    if(allocated(monoPREC%DMUMPS%offset_list))then
      deallocate(monoPREC%DMUMPS%offset_list)
    endif

    if(allocated(monoPREC%DMUMPS%offset_counts))then
      deallocate(monoPREC%DMUMPS%offset_counts)
    endif

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
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        jn = monoMAT%CSR%item(j)
        if(is_self .and. monoMAT%N < jn) cycle aa
        NZ = NZ + 1
      enddo
    enddo aa
  end subroutine monolis_precond_mumps_get_nz

  subroutine monolis_precond_mumps_get_loc(monoMAT, monoPREC, monoCOM, IRN, JCN, A, is_self)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoPREC
    type(monolis_com) :: monoCOM
    integer(kint) :: i, j, in, jS, jE
    integer(kint) :: idof, jdof, jn, kn, NDOF, nprocs, myrank, i2, jn2, offset
    integer(kint), allocatable :: offset_id(:), offset_rank(:)
    integer(kint), pointer :: IRN(:), JCN(:)
    real(kdouble) :: A(:)
    logical :: is_self

    NDOF = monoMAT%NDOF
    nprocs = monolis_mpi_get_local_comm_size(monoCOM%comm)
    myrank = monolis_mpi_get_local_my_rank(monoCOM%comm)

    offset = monoPREC%DMUMPS%offset_list(monolis_mpi_get_local_my_rank(monoCOM%comm)+1)
    allocate(offset_id(monoMAT%NP), source = 0)
    allocate(offset_rank(monoMAT%NP), source = -1)

    do i = 1, monoMAT%NP
      offset_id(i) = i
    enddo
    call monolis_mpi_update_I(monoCOM, 1, offset_id)

    do i = 1, monoCOM%recv_n_neib
      jS = monoCOM%recv_index(i) + 1
      jE = monoCOM%recv_index(i + 1)
      do j = jS, jE
        in = monoCOM%recv_item(j)
        offset_rank(in) = monoCOM%recv_neib_pe(i)
      enddo
    enddo

    do i = 1, monoMAT%NP
      if(offset_rank(i) == -1)then
        offset_id(i) = offset_id(i) + monoPREC%DMUMPS%offset_list(myrank+1)/monoMAT%NDOF
      else
        offset_id(i) = offset_id(i) + monoPREC%DMUMPS%offset_list(offset_rank(i)+1)/monoMAT%NDOF
      endif
    enddo

    in = 0
    aa:do i = 1, monoMAT%N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        jn = monoMAT%CSR%item(j)
        if(is_self .and. monoMAT%N < jn) cycle aa
        do idof = 1, NDOF
        do jdof = 1, NDOF
          kn = NDOF*NDOF*(j-1) + NDOF*(idof-1) + jdof
          in = in + 1
          A(in) = monoMAT%R%A(kn)
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

  subroutine monolis_precond_mumps_set_rhs(N, X, RHS, comm, &
    & offset_list, offset_counts, is_self)
    implicit none
    integer(kint) :: N, comm
    real(kdouble) :: X(:), RHS(:)
    integer(kint) :: offset_list(:)
    integer(kint) :: offset_counts(:)
    logical :: is_self

    if(monolis_mpi_get_local_comm_size(comm) == 1 .or. is_self)then
      RHS = X(1:N)
    else
      call monolis_gatherv_R(X, N, &
        RHS, offset_counts, offset_list, &
        0, comm)
    endif
  end subroutine monolis_precond_mumps_set_rhs

  subroutine monolis_precond_mumps_get_rhs(N, RHS, Y, comm, &
    & offset_list, offset_counts, is_self)
    implicit none
    integer(kint) :: N, comm
    real(kdouble) :: Y(:), RHS(:)
    integer(kint) :: offset_list(:)
    integer(kint) :: offset_counts(:)
    logical :: is_self

    if(monolis_mpi_get_local_comm_size(comm) == 1 .or. is_self)then
      Y(1:N) = RHS
    else
      call monolis_scatterv_R( &
        RHS, offset_counts, offset_list, &
        Y, N, &
        0, comm)
    endif
  end subroutine monolis_precond_mumps_get_rhs

  subroutine monolis_precond_mumps_get_offset(N_loc, comm, offset_list, offset_counts)
    implicit none
    integer(kint) :: nprocs, myrank, i, N_loc, comm
    integer(kint), allocatable :: offset_list(:)
    integer(kint), allocatable :: offset_counts(:)

    nprocs = monolis_mpi_get_local_comm_size(comm)
    myrank = monolis_mpi_get_local_my_rank(comm)

    allocate(offset_list(nprocs), source = 0)
    allocate(offset_counts(nprocs), source = 0)

    if(1 < nprocs)then
      call monolis_allgather_I1(N_loc, offset_counts, comm)
    endif

    offset_list(1) = 0
    do i = 1, nprocs - 1
      offset_list(i+1) = offset_list(i) + offset_counts(i)
    enddo
  end subroutine monolis_precond_mumps_get_offset
end module mod_monolis_precond_mumps
