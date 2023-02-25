!> ソルバパラメータの定義
module mod_monolis_def_solver
  use mod_monolis_utils
  implicit none

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
  integer(kint), parameter :: monolis_prec_LU     = 8
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
  & "LU    ", &
  & "MUMPS ", &
  & "ROM   ", &
  & "MF    ", &
  & "MUMPSL"/)

  integer(kint), parameter :: monolis_prm_Iarray_size = 100
  integer(kint), parameter :: monolis_prm_Rarray_size = 100
  integer(kint), parameter :: monolis_prm_Iarray(monolis_prm_Iarray_size)
  real(kdouble), parameter :: monolis_prm_Rarray(monolis_prm_Rarray_size)

  integer(kint), parameter :: monolis_prm_method = 1
  integer(kint), parameter :: monolis_prm_precond = 2
  integer(kint), parameter :: monolis_prm_max_iter = 3
  integer(kint), parameter :: monolis_prm_cur_iter = 4
  integer(kint), parameter :: monolis_prm_ierr = 5
  integer(kint), parameter :: monolis_prm_is_scaling = 6
  integer(kint), parameter :: monolis_prm_is_reordering = 7
  integer(kint), parameter :: monolis_prm_is_init_x = 8
  integer(kint), parameter :: monolis_prm_is_sym_matrix = 9
  integer(kint), parameter :: monolis_prm_is_debug = 10
  integer(kint), parameter :: monolis_prm_is_measurement = 11
  integer(kint), parameter :: monolis_prm_is_check_diag = 12
  integer(kint), parameter :: monolis_prm_is_prec_stored = 13
  integer(kint), parameter :: monolis_prm_show_iterlog = 14
  integer(kint), parameter :: monolis_prm_show_time = 15
  integer(kint), parameter :: monolis_prm_show_time_statistics = 16
  integer(kint), parameter :: monolis_prm_show_summary = 17

  integer(kint), parameter :: monolis_prm_tol = 1
  integer(kint), parameter :: monolis_prm_cur_resid = 2
  integer(kint), parameter :: monolis_time_sol = 3
  integer(kint), parameter :: monolis_time_prep = 4
  integer(kint), parameter :: monolis_time_spmv = 5
  integer(kint), parameter :: monolis_time_dotp = 6
  integer(kint), parameter :: monolis_time_prec = 7
  integer(kint), parameter :: monolis_time_comm_dotp = 8
  integer(kint), parameter :: monolis_time_comm_spmv = 9

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
  end type monolis_prm

contains

  subroutine monolis_prm_initialize(monoPRM, fname_in)
    implicit none
    type(monolis_prm) :: monoPRM
    character(*) :: fname_in

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
end module mod_monolis_def_solver