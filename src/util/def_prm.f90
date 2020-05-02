module mod_monolis_prm
  implicit none

  integer(4), parameter :: kint    = 4
  integer(4), parameter :: kdouble = 8

  integer(4), parameter :: monolis_iter_CG       = 1
  integer(4), parameter :: monolis_iter_GropCG   = 2
  integer(4), parameter :: monolis_iter_PipeCG   = 3
  integer(4), parameter :: monolis_iter_PipeCR   = 4
  integer(4), parameter :: monolis_iter_BiCGSTAB = 5
  integer(4), parameter :: monolis_iter_PipeBiCGSTAB = 6
  integer(4), parameter :: monolis_iter_BiCGSTAB_noprec = 7
  integer(4), parameter :: monolis_iter_CABiCGSTAB_noprec = 8
  integer(4), parameter :: monolis_iter_PipeBiCGSTAB_noprec = 9
  integer(4), parameter :: monolis_iter_SOR      = 10
  integer(4), parameter :: monolis_iter_IR       = 11

  integer(4), parameter :: monolis_prec_NONE   = 0
  integer(4), parameter :: monolis_prec_DIAG   = 1
  integer(4), parameter :: monolis_prec_ILU    = 2
  integer(4), parameter :: monolis_prec_JACOBI = 3
  integer(4), parameter :: monolis_prec_SOR    = 4
  integer(4), parameter :: monolis_prec_SAINV  = 5
  integer(4), parameter :: monolis_prec_RIF    = 6
  integer(4), parameter :: monolis_prec_SPIKE  = 7
  integer(4), parameter :: monolis_prec_DIRECT = 8
  integer(4), parameter :: monolis_prec_MUMPS  = 9

  !enum, bind(c)
  !  enumerator :: BoundaryType_Dirichlet = 1
  !  enumerator :: BoundaryType_Neumann
  !  enumerator :: BoundaryType_Robin
  !end enum

  character*24, dimension(11) :: monolis_str_iter = (/&
  & "CG                 ", &
  & "GropCG             ", &
  & "PipeCG             ", &
  & "PipeCR             ", &
  & "BiCGSTAB           ", &
  & "PipeBiCGSTAB       ", &
  & "BiCGSTAB_noprec    ", &
  & "CABiCGSTAB_noprec  ", &
  & "PipeBiCGSTAB_noprec", &
  & "SOR                ", &
  & "IR                 "/)
  character*24, dimension(0:9)  :: monolis_str_prec = (/&
  & "None  ", &
  & "Diag  ", &
  & "ILU   ", &
  & "Jacobi", &
  & "SOR   ", &
  & "SAINV ", &
  & "RIF   ", &
  & "SPIKE ", &
  & "Direct", &
  & "MUMPS "/)

  type monolis_prm
    integer(kind=kint) :: method = 1
    integer(kind=kint) :: precond = 1
    integer(kind=kint) :: maxiter = 1000
    integer(kind=kint) :: curiter = 0
    integer(kind=kint) :: ierr = -1
    real(kind=kdouble) :: tol = 1.0d-8
    real(kind=kdouble) :: curresid = 0.0d0
    logical :: is_scaling    = .true.
    logical :: is_reordering = .true.
    logical :: is_init_x     = .true.
    logical :: is_sym_matrix = .true.
    logical :: is_debug      = .false.
    logical :: is_check_diag = .true.
    logical :: show_iterlog  = .true.
    logical :: show_time     = .true.
    logical :: show_summary  = .true.
    !> time: tsol = tspmv + tprec + tcomm + others
    real(kind=kdouble) :: tsol  = 0.0d0
    real(kind=kdouble) :: tprep = 0.0d0
    real(kind=kdouble) :: tspmv = 0.0d0
    real(kind=kdouble) :: tdotp = 0.0d0
    real(kind=kdouble) :: tprec = 0.0d0
    real(kind=kdouble) :: tcomm = 0.0d0
  end type monolis_prm

contains

  subroutine monolis_prm_initialize(monoPRM)
    implicit none
    type(monolis_prm) :: monoPRM

    monoPRM%method = 1
    monoPRM%precond = 1
    monoPRM%curiter = 0
    monoPRM%maxiter = 1000
    monoPRM%ierr = -1
    monoPRM%tol = 1.0d-8
    monoPRM%curresid = 0.0d0
    monoPRM%is_scaling    = .true.
    monoPRM%is_reordering = .true.
    monoPRM%is_init_x     = .true.
    monoPRM%is_debug      = .false.
    monoPRM%show_iterlog  = .true.
    monoPRM%show_time     = .true.
    monoPRM%show_summary  = .true.
  end subroutine monolis_prm_initialize

  subroutine monolis_prm_finalize(monoPRM)
    implicit none
    type(monolis_prm) :: monoPRM

  end subroutine monolis_prm_finalize

end module mod_monolis_prm