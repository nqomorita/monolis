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
  public :: monolis_finalize
  public :: monolis_timer_initialize
  public :: monolis_timer_finalize
  public :: monolis_check_diagonal
  public :: monolis_set_debug
  public :: monolis_debug_header
  public :: monolis_debug_int
  public :: monolis_debug_char
  public :: monolis_debug_logical
  public :: monolis_warning_header
  public :: monolis_get_time
  public :: monolis_get_time_sync
  public :: monolis_error_stop

  integer(kint), save :: myrank = 0
  integer(kint), save :: mycomm
  logical, save :: is_debug

contains

  subroutine monolis_global_initialize()
    implicit none
    call monolis_mpi_initialize()
    myrank = monolis_global_myrank()
    mycomm = monolis_global_comm()
  end subroutine monolis_global_initialize

  subroutine monolis_global_finalize()
    implicit none
    call monolis_mpi_finalize()
  end subroutine monolis_global_finalize

  subroutine monolis_initialize(monolis)
    implicit none
    type(monolis_structure) :: monolis

    call monolis_prm_initialize(monolis%PRM)
    call monolis_com_initialize(monolis%COM)
    call monolis_mat_initialize(monolis%MAT)
    call monolis_com_input_comm_table(monolis%COM)
  end subroutine monolis_initialize

  subroutine monolis_finalize(monolis)
    implicit none
    type(monolis_structure) :: monolis

    call monolis_prm_finalize(monolis%PRM)
    call monolis_com_finalize(monolis%COM)
    call monolis_mat_finalize(monolis%MAT)
  end subroutine monolis_finalize

  subroutine monolis_timer_initialize(monoPRM, monoCOM)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM

    is_debug = monoPRM%is_debug
    call monolis_debug_header("monolis_timer_initialize")

    call monolis_barrier(monoCOM)
    monoPRM%tsol  = monolis_get_time()
    monoPRM%tprep = 0.0d0
    monoPRM%tspmv = 0.0d0
    monoPRM%tdotp = 0.0d0
    monoPRM%tprec = 0.0d0
    monoPRM%tcomm_dotp = 0.0d0
    monoPRM%tcomm_spmv = 0.0d0
  end subroutine monolis_timer_initialize

  subroutine monolis_timer_finalize(monoPRM, monoCOM)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    real(kdouble) :: t1
    logical :: is_output

    call monolis_debug_header("monolis_timer_finalize")

    t1 = monolis_get_time()
    monoPRM%tsol = t1 - monoPRM%tsol

    if(monoPRM%show_summary .and. monoCOM%myrank == 0)then
      write(*,"(a,i10)")" ** monolis converge iter:", monoPRM%curiter
      write(*,"(a,1p4e10.3)")" ** monolis rel. residual:", monoPRM%curresid
    endif

    is_output = monoPRM%show_summary .or. monoPRM%show_time
    if(is_output .and. monoCOM%myrank == 0)then
      write(*,"(a,1p4e10.3)")" ** monolis solution time:", monoPRM%tsol
    endif

    if(monoPRM%show_time .and. monoCOM%myrank == 0)then
      write(*,"(a,1p4e10.3)")"  - solution/prepost time:", monoPRM%tprep
      write(*,"(a,1p4e10.3)")"  - solution/SpMV    time:", monoPRM%tspmv
      write(*,"(a,1p4e10.3)")"  - solution/inner p time:", monoPRM%tdotp
      write(*,"(a,1p4e10.3)")"  - solution/precond time:", monoPRM%tprec
      write(*,"(a,1p4e10.3)")"  - (solution/comm dotp) :", monoPRM%tcomm_dotp
      write(*,"(a,1p4e10.3)")"  - (solution/comm spmv) :", monoPRM%tcomm_spmv
    endif
  end subroutine monolis_timer_finalize

  function monolis_get_input_filename(input)
    implicit none
    integer(kint) :: comm_size, myrank, pos
    character(*) :: input
    character :: cnum*6
    character :: monolis_get_input_filename*128

    comm_size = monolis_global_commsize()
    monolis_get_input_filename = trim(input)

    if(comm_size > 1)then
      myrank = monolis_global_myrank()
      write(cnum,"(i0)") myrank

      !> for node.dat
      pos = index(monolis_get_input_filename, "node.dat", back = .true.)
      if(pos > 0)then
        monolis_get_input_filename(pos:pos+16) = "parted/node.dat."
        monolis_get_input_filename = trim(monolis_get_input_filename) // trim(cnum)
      endif

      !> for elem.dat
      pos = index(monolis_get_input_filename, "elem.dat", back = .true.)
      if(pos > 0)then
        monolis_get_input_filename(pos:pos+16) = "parted/elem.dat."
        monolis_get_input_filename = trim(monolis_get_input_filename) // trim(cnum)
      endif
    endif
  end function monolis_get_input_filename

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

  function monolis_get_time_sync()
    implicit none
    real(kdouble) :: monolis_get_time_sync, t1

#ifdef WITH_MPI
    call monolis_barrier(mycomm)
    monolis_get_time_sync = MPI_Wtime()
#else
    call cpu_time(t1)
    monolis_get_time_sync = t1
#endif
  end function monolis_get_time_sync

  subroutine monolis_check_diagonal(monoPRM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, k, jS, jE, in, kn, NP, NDOF, NDOF2
    real(kdouble) :: t1, t2

    if(.not. monoPRM%is_check_diag) return
    call monolis_debug_header("monolis_check_diagonal")
    t1 = monolis_get_time()

    NP =  monoMAT%NP
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF

    do i = 1, NP
      jS = monoMAT%index(i-1) + 1
      jE = monoMAT%index(i)
      do j = jS, jE
        in = monoMAT%item(j)
        if(i == in)then
          do k = 1, NDOF
            kn = NDOF2*(j-1) + (NDOF+1)*(k-1) + 1
            if(monoMAT%A(kn) == 0.0d0)then
              if(myrank == 0) write(*,"(a,i8,a,i8)")" ** monolis error: zero diagonal at node:", i, " , dof: ", k
              stop
            endif
          enddo
        endif
      enddo
    enddo

    t2 = monolis_get_time()
    monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_check_diagonal

  !> set parameter section
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

  subroutine monolis_param_set_is_scaling(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_scaling = param
  end subroutine monolis_param_set_is_scaling

  subroutine monolis_param_set_is_reordering(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_reordering = param
  end subroutine monolis_param_set_is_reordering

  subroutine monolis_param_set_is_init_x(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_init_x = param
  end subroutine monolis_param_set_is_init_x

  subroutine monolis_param_set_is_sym_matrix(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_sym_matrix = param
  end subroutine monolis_param_set_is_sym_matrix

  subroutine monolis_param_set_is_debug(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_debug = param
  end subroutine monolis_param_set_is_debug

  subroutine monolis_param_set_is_check_diag(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    logical :: param
    monolis%PRM%is_check_diag = param
  end subroutine monolis_param_set_is_check_diag

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

  !> debug section
  subroutine monolis_set_debug(flag)
    implicit none
    logical :: flag
    is_debug = flag
  end subroutine monolis_set_debug

  subroutine monolis_debug_header(header)
    implicit none
    character(*) :: header

    if(.not. is_debug) return
    if(myrank == 0) write(*,"(a)")"** monolis debug: "//trim(header)
  end subroutine monolis_debug_header

  subroutine monolis_debug_int(header, n)
    implicit none
    integer(kint) :: n
    character(*) :: header

    if(.not. is_debug) return
    if(myrank == 0) write(*,"(a,i12)")"** monolis debug: "//trim(header)//": ", n
  end subroutine monolis_debug_int

  subroutine monolis_debug_char(header, char)
    implicit none
    character(*) :: header, char

    if(.not. is_debug) return
    if(myrank == 0) write(*,"(a,a)")"** monolis debug: "//trim(header)//": ", trim(char)
  end subroutine monolis_debug_char

  subroutine monolis_debug_logical(header, l)
    implicit none
    logical :: l
    character(*) :: header

    if(.not. is_debug) return
    if(myrank == 0) write(*,"(a,l)")"** monolis debug: "//trim(header)//": ", l
  end subroutine monolis_debug_logical

  subroutine monolis_warning_header(header)
    implicit none
    character(*) :: header

    if(myrank == 0) write(*,"(a)")"** monolis warning: "//trim(header)
  end subroutine monolis_warning_header

  subroutine monolis_error_stop()
    implicit none
    error stop monolis_fail
  end subroutine monolis_error_stop
end module mod_monolis_util