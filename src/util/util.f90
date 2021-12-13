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
  public :: monolis_timer_initialize
  !public :: monolis_timer_finalize
  public :: monolis_get_input_filename
  public :: monolis_check_diagonal
  public :: monolis_set_debug
  public :: monolis_debug_header
  public :: monolis_debug_int
  public :: monolis_debug_real
  public :: monolis_debug_char
  public :: monolis_debug_logical
  public :: monolis_warning_header
  public :: monolis_get_time
  public :: monolis_get_time_sync
  public :: monolis_error_stop
  public :: monolis_get_penalty_value
  public :: monolis_param_set_method
  public :: monolis_param_set_precond
  public :: monolis_param_set_maxiter
  public :: monolis_param_set_tol
  public :: monolis_param_set_is_scaling
  public :: monolis_param_set_is_reordering
  public :: monolis_param_set_is_init_x
  public :: monolis_param_set_is_sym_matrix
  public :: monolis_param_set_is_debug
  public :: monolis_param_set_is_check_diag
  public :: monolis_param_set_show_iterlog
  public :: monolis_param_set_show_time
  public :: monolis_param_set_show_summary
  public :: monolis_get_internal_elem_1d_bool

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

  subroutine monolis_timer_initialize(monoPRM, monoCOM)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM

    is_debug = monoPRM%is_debug
    call monolis_debug_header("monolis_timer_initialize")

    call monolis_barrier(monoCOM%comm)
    monoPRM%tsol  = monolis_get_time()
    monoPRM%tprep = 0.0d0
    monoPRM%tspmv = 0.0d0
    monoPRM%tdotp = 0.0d0
    monoPRM%tprec = 0.0d0
    monoPRM%tcomm_dotp = 0.0d0
    monoPRM%tcomm_spmv = 0.0d0
  end subroutine monolis_timer_initialize

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
      output_dir = "parted/"
      if(present(fname_dir)) output_dir = trim(fname_dir)//"/parted/"
      write(cnum,"(i0)") myrank
      monolis_get_input_filename = trim(output_dir)//trim(fname_in)//"."//trim(cnum)
    else
      monolis_get_input_filename = trim(fname_in)
      if(present(fname_dir)) monolis_get_input_filename = trim(fname_dir)//"/"//trim(fname_in)
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

  function monolis_get_penalty_value(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, k, jS, jE, in, kn, NP, NDOF, NDOF2
    real(kdouble) :: monolis_get_penalty_value, max

    call monolis_debug_header("monolis_get_penalty_value")

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
    call monolis_debug_header("monolis_check_diagonal")
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

  subroutine monolis_param_set_tol(monolis, param)
    implicit none
    type(monolis_structure) :: monolis
    real(kdouble) :: param
    monolis%PRM%tol = param
  end subroutine monolis_param_set_tol

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

  !> debug section
  subroutine monolis_print_char(header, char)
    implicit none
    character(*) :: header, char
    if(myrank == 0) write(*,"(a,a)")">> "//trim(header)//": ", trim(char)
  end subroutine monolis_print_char

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

  subroutine monolis_debug_real(header, v)
    implicit none
    real(kdouble) :: v
    character(*) :: header

    if(.not. is_debug) return
    if(myrank == 0) write(*,"(a,1pe12.4)")"** monolis debug: "//trim(header)//": ", v
  end subroutine monolis_debug_real

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
