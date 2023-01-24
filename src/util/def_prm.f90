module mod_monolis_prm
  implicit none

  integer(4), parameter :: kint    = 4
  integer(4), parameter :: kdouble = 8
  integer(kint), parameter :: monolis_success = 0
  integer(kint), parameter :: monolis_fail = 1
  integer(kint), parameter :: monolis_charlen = 1024

  integer(kint), parameter :: monolis_iter_CG       = 1
  integer(kint), parameter :: monolis_iter_GropCG   = 2
  integer(kint), parameter :: monolis_iter_PipeCG   = 3
  integer(kint), parameter :: monolis_iter_PipeCR   = 4
  integer(kint), parameter :: monolis_iter_BiCGSTAB = 5
  integer(kint), parameter :: monolis_iter_PipeBiCGSTAB = 6
  integer(kint), parameter :: monolis_iter_BiCGSTAB_noprec = 7
  integer(kint), parameter :: monolis_iter_CABiCGSTAB_noprec = 8
  integer(kint), parameter :: monolis_iter_PipeBiCGSTAB_noprec = 9
  integer(kint), parameter :: monolis_iter_SOR      = 10
  integer(kint), parameter :: monolis_iter_IR       = 11
  integer(kint), parameter :: monolis_iter_GMRES    = 12

  integer(kint), parameter :: monolis_prec_NONE   = 0
  integer(kint), parameter :: monolis_prec_DIAG   = 1
  integer(kint), parameter :: monolis_prec_ILU    = 2
  integer(kint), parameter :: monolis_prec_JACOBI = 3
  integer(kint), parameter :: monolis_prec_SOR    = 4
  integer(kint), parameter :: monolis_prec_SAINV  = 5
  integer(kint), parameter :: monolis_prec_RIF    = 6
  integer(kint), parameter :: monolis_prec_SPIKE  = 7
  integer(kint), parameter :: monolis_prec_DIRECT = 8
  integer(kint), parameter :: monolis_prec_MUMPS  = 9
  integer(kint), parameter :: monolis_prec_ROM    = 10
  integer(kint), parameter :: monolis_prec_MF     = 11
  integer(kint), parameter :: monolis_prec_MUMPS_LOCAL = 12

  character*24, dimension(12) :: monolis_str_iter = (/&
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
  & "IR                 ", &
  & "GMRES              "/)

  character*24, dimension(0:12)  :: monolis_str_prec = (/&
  & "None  ", &
  & "Diag  ", &
  & "ILU   ", &
  & "Jacobi", &
  & "SOR   ", &
  & "SAINV ", &
  & "RIF   ", &
  & "SPIKE ", &
  & "Direct", &
  & "MUMPS ", &
  & "ROM   ", &
  & "MF    ", &
  & "MUMPS "/)

  type monolis_prm
    integer(kint) :: method = 1
    integer(kint) :: precond = 1
    integer(kint) :: maxiter = 1000
    integer(kint) :: curiter = 0
    integer(kint) :: ierr = -1
    real(kdouble) :: tol = 1.0d-8
    real(kdouble) :: curresid = 0.0d0
    logical :: is_scaling    = .false.
    logical :: is_reordering = .false.
    logical :: is_init_x     = .true.
    logical :: is_sym_matrix = .true.
    logical :: is_debug      = .false.
    logical :: is_measurement= .false.
    logical :: is_check_diag = .false.
    logical :: is_prec_stored= .false.
    logical :: show_iterlog  = .true.
    logical :: show_time     = .true.
    logical :: show_time_statistics = .false.
    logical :: show_summary  = .true.
    !> time: tsol
    real(kdouble) :: tsol  = 0.0d0
    real(kdouble) :: tprep = 0.0d0
    real(kdouble) :: tspmv = 0.0d0
    real(kdouble) :: tdotp = 0.0d0
    real(kdouble) :: tprec = 0.0d0
    real(kdouble) :: tcomm_dotp = 0.0d0
    real(kdouble) :: tcomm_spmv = 0.0d0
    !> input dir
    !> character(monolis_charlen) :: input_file_dir
  end type monolis_prm

contains

  subroutine monolis_prm_initialize(monoPRM, fname_in)
    implicit none
    type(monolis_prm) :: monoPRM
    character(*) :: fname_in

    !> monoPRM%input_file_dir = trim(fname_in)
    monoPRM%method = 1
    monoPRM%precond = 1
    monoPRM%curiter = 0
    monoPRM%maxiter = 1000
    monoPRM%ierr = -1
    monoPRM%tol = 1.0d-8
    monoPRM%curresid = 0.0d0
    monoPRM%is_scaling    = .false.
    monoPRM%is_reordering = .false.
    monoPRM%is_init_x     = .true.
    monoPRM%is_sym_matrix = .true.
    monoPRM%is_debug      = .false.
    monoPRM%is_measurement= .false.
    monoPRM%is_check_diag = .true.
    monoPRM%show_iterlog  = .true.
    monoPRM%show_time     = .true.
    monoPRM%show_summary  = .true.
  end subroutine monolis_prm_initialize

  subroutine monolis_prm_finalize(monoPRM)
    implicit none
    type(monolis_prm) :: monoPRM

  end subroutine monolis_prm_finalize

end module mod_monolis_prm