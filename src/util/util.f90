module mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

  type monolis_structure
    type(monolis_prm) :: PRM
    type(monolis_com) :: COM
    type(monolis_mat) :: MAT
  end type monolis_structure

  private
  public :: monolis_structure
  public :: monolis_global_initialize
  public :: monolis_global_finalize
  public :: monolis_initialize
  public :: monolis_initialize_entire
  public :: monolis_finalize
  public :: monolis_get_time
  public :: monolis_get_time_sync
  public :: monolis_get_input_filename
  public :: monolis_check_diagonal
  public :: monolis_get_penalty_value
  public :: monolis_get_internal_elem_1d_bool

contains

  subroutine monolis_global_initialize()
    implicit none
    call monolis_mpi_initialize()
  end subroutine monolis_global_initialize

  subroutine monolis_global_finalize()
    implicit none
    call monolis_mpi_finalize()
  end subroutine monolis_global_finalize

  subroutine monolis_initialize(monolis, fname_in)
    implicit none
    type(monolis_structure) :: monolis
    character(*) :: fname_in

    call monolis_prm_initialize(monolis%PRM, fname_in)
    call monolis_com_initialize(monolis%COM)
    call monolis_mat_initialize(monolis%MAT)
    call monolis_com_input_comm_table(monolis%COM, fname_in)
  end subroutine monolis_initialize

  subroutine monolis_initialize_entire(monolis, fname_in)
    implicit none
    type(monolis_structure) :: monolis
    character(*) :: fname_in

    call monolis_prm_initialize(monolis%PRM, fname_in)
    call monolis_com_initialize(monolis%COM, .true.)
    call monolis_mat_initialize(monolis%MAT)
    call monolis_com_input_comm_table(monolis%COM, fname_in)
  end subroutine monolis_initialize_entire

  subroutine monolis_finalize(monolis)
    implicit none
    type(monolis_structure) :: monolis

    call monolis_prm_finalize(monolis%PRM)
    call monolis_com_finalize(monolis%COM)
    call monolis_mat_finalize(monolis%MAT)
  end subroutine monolis_finalize

  function monolis_get_time()
    implicit none
    real(kdouble) :: monolis_get_time, t1

#ifdef WITH_MPI
    monolis_get_time = MPI_Wtime()
#else
    call cpu_time(t1)
    monolis_get_time = t1
#endif
  end function monolis_get_time

  function monolis_get_time_sync(comm)
    implicit none
    real(kdouble) :: monolis_get_time_sync, t1
    integer(kint) :: mycomm
    integer(kint), optional :: comm

#ifdef WITH_MPI
    if(.not. present(comm))then
      mycomm = monolis_global_comm()
    endif

    call monolis_barrier_(mycomm)
    monolis_get_time_sync = MPI_Wtime()
#else
    call cpu_time(t1)
    monolis_get_time_sync = t1
#endif
  end function monolis_get_time_sync

  function monolis_get_input_filename(fname_in, fname_dir)
    implicit none
    integer(kint) :: comm_size, myrank
    character(*) :: fname_in
    character :: cnum*6, output_dir*100
    character(monolis_charlen) :: monolis_get_input_filename
    character(*), optional :: fname_dir

    comm_size = monolis_global_commsize()
    myrank = monolis_global_myrank()
    if(comm_size > 1)then
      output_dir = "parted.0/"
      if(present(fname_dir)) output_dir = trim(fname_dir)//"/parted.0/"
      write(cnum,"(i0)") myrank
      monolis_get_input_filename = trim(output_dir)//trim(fname_in)//"."//trim(cnum)
    else
      monolis_get_input_filename = trim(fname_in)
      if(present(fname_dir)) monolis_get_input_filename = trim(fname_dir)//"/"//trim(fname_in)
    endif
  end function monolis_get_input_filename

  function monolis_get_penalty_value(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, k, jS, jE, in, kn, NP, NDOF, NDOF2
    real(kdouble) :: monolis_get_penalty_value, max

    NP =  monoMAT%NP
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    max = 0.0d0

    do i = 1, NP
      jS = monoMAT%index(i-1) + 1
      jE = monoMAT%index(i)
      do j = jS, jE
        in = monoMAT%item(j)
        if(i == in)then
          do k = 1, NDOF
            kn = NDOF2*(j-1) + (NDOF+1)*(k-1) + 1
            if(max < monoMAT%A(kn)) max = monoMAT%A(kn)
          enddo
        endif
      enddo
    enddo
    monolis_get_penalty_value = max
  end function monolis_get_penalty_value

  subroutine monolis_check_diagonal(monoPRM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, k, jS, jE, in, kn, N, NDOF, NDOF2
    real(kdouble) :: t1, t2

    if(.not. monoPRM%is_check_diag) return
    t1 = monolis_get_time()

    N =  monoMAT%N
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF

    do i = 1, N
      jS = monoMAT%index(i-1) + 1
      jE = monoMAT%index(i)
      do j = jS, jE
        in = monoMAT%item(j)
        if(i == in)then
          do k = 1, NDOF
            kn = NDOF2*(j-1) + (NDOF+1)*(k-1) + 1
            if(monoMAT%A(kn) == 0.0d0)then
              write(*,"(a,i8,a,i8)")" ** monolis error: zero diagonal at node:", i, " , dof: ", k
              stop
            endif
          enddo
        endif
      enddo
    enddo

    t2 = monolis_get_time()
    monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_check_diagonal

  subroutine monolis_get_internal_elem_1d_bool(monolis, nelem, list)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: nelem, i
    logical :: list(:)

    if(monolis_global_commsize() == 1)then
      list = .true.
    else
      list = .false.
      do i = 1, monolis%COM%internal_nelem
        list(i) = .true.
      enddo
    endif
  end subroutine monolis_get_internal_elem_1d_bool

end module mod_monolis_util
