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
  public :: monolis_initialize
  public :: monolis_finalize
  public :: monolis_timer_initialize
  public :: monolis_timer_finalize
  public :: monolis_check_diagonal
  public :: monolis_debug_header
  public :: monolis_get_time

  integer(kind=kint), save :: myrank

contains

  subroutine monolis_initialize(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_prm_initialize(monoPRM)
    call monolis_com_initialize(monoCOM)
    call monolis_mat_initialize(monoMAT)
    myrank = monoCOM%myrank
  end subroutine monolis_initialize

  subroutine monolis_finalize(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_prm_finalize(monoPRM)
    call monolis_com_finalize(monoCOM)
    call monolis_mat_finalize(monoMAT)
  end subroutine monolis_finalize

  subroutine monolis_timer_initialize(monoPRM)
    implicit none
    type(monolis_prm) :: monoPRM

    if(monoPRM%is_debug) call monolis_debug_header("monolis_timer_initialize")

    monoPRM%tsol  = monolis_get_time()
    monoPRM%tprep = 0.0d0
    monoPRM%tspmv = 0.0d0
    monoPRM%tdotp = 0.0d0
    monoPRM%tprec = 0.0d0
    monoPRM%tcomm = 0.0d0
  end subroutine monolis_timer_initialize

  subroutine monolis_timer_finalize(monoPRM, monoCOM)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    real(kind=kdouble) :: t1
    logical :: is_output

    if(monoPRM%is_debug) call monolis_debug_header("monolis_timer_finalize")

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
      !write(*,"(a,1p4e10.3)")"  - solution/comm    time:", monoPRM%tcomm
    endif
  end subroutine monolis_timer_finalize

  function monolis_get_time()
    implicit none
    real(kind=kdouble) :: monolis_get_time, t1

#ifdef WITH_MPI
    monolis_get_time = MPI_Wtime()
#else
    call cpu_time(t1)
    monolis_get_time = t1
#endif
  end function monolis_get_time

  subroutine monolis_check_diagonal(monoPRM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, k, jS, jE, in, kn, NP, NDOF, NDOF2
    real(kind=kdouble) :: t1, t2

    if(.not. monoPRM%is_check_diag) return
    if(monoPRM%is_debug) call monolis_debug_header("monolis_check_diagonal")
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

  subroutine monolis_debug_header(header)
    implicit none
    character(*) :: header

    if(myrank == 0) write(*,"(a)")"** monolis debug: "//trim(header)
  end subroutine monolis_debug_header

  subroutine monolis_debug_equal_R(header, a, b)
    implicit none
    character(*) :: header
    real(kind=kdouble) :: a, b

    if(a /= b)then
      if(myrank == 0) write(*,"(a,1pe12.5,a,1pe12.5,a)")"** monolis debug: ", a, "is not equal", b, " at "//trim(header)
    endif
  end subroutine monolis_debug_equal_R

end module mod_monolis_util