!> 行列構造体の定義
module mod_monolis_def_mat
  use mod_monolis_utils
  implicit none

  !> 行列構造体（実数型）
  type monolis_mat_val_R
    !> 全行列値
    real(kdouble), allocatable :: A(:)
    !> 行列値（狭義上三角）
    real(kdouble), allocatable :: U(:)
    !> 行列値（対角成分）
    real(kdouble), allocatable :: D(:)
    !> 行列値（狭義下三角）
    real(kdouble), allocatable :: L(:)
    !> 解ベクトル
    real(kdouble), allocatable :: X(:)
    !> 右辺ベクトル
    real(kdouble), allocatable :: B(:)
  end type monolis_mat_val_R

  !> 行列構造体（複素数型）
  type monolis_mat_val_C
    !> 全行列値
    complex(kdouble), allocatable :: A(:)
    !> 行列値（狭義上三角）
    complex(kdouble), allocatable :: U(:)
    !> 行列値（対角成分）
    complex(kdouble), allocatable :: D(:)
    !> 行列値（狭義下三角）
    complex(kdouble), allocatable :: L(:)
    !> 解ベクトル
    complex(kdouble), allocatable :: X(:)
    !> 右辺ベクトル
    complex(kdouble), allocatable :: B(:)
  end type monolis_mat_val_C

  !> 行列構造体（セパレート CSR 構造）
  type monolis_mat_separated_CSR
    !> index 配列（狭義上三角）
    integer(kint), allocatable :: indexU(:)
    !> item 配列（狭義上三角）
    integer(kint), allocatable :: itemU(:)
    !> index 配列（狭義下三角）
    integer(kint), allocatable :: indexL(:)
    !> item 配列（狭義下三角）
    integer(kint), allocatable :: itemL(:)
  end type monolis_mat_separated_CSR

  !> 行列構造体（CSR 構造）
  type monolis_mat_CSR
    !> index 配列
    integer(kint), allocatable :: index(:)
    !> item 配列
    integer(kint), allocatable :: item(:)
  end type monolis_mat_CSR

  !> 行列構造体（CSC 構造）
  type monolis_mat_CSC
    !> index 配列
    integer(kint), allocatable :: index(:)
    !> item 配列
    integer(kint), allocatable :: item(:)
    !> CSR 形式に対する行列値の置換ベクトル
    integer(kint), allocatable :: perm(:)
  end type monolis_mat_CSC

  !> 行列構造体（セパレート CSR 構造）
  !type monolis_mat_lattice
  !end type monolis_mat_lattice

  !> 行列構造体
  type monolis_mat
    !> 内部自由度数
    integer(kint) :: N
    !> 全自由度数
    integer(kint) :: NP
    !> 1 ブロックの自由度
    integer(kint) :: NDOF
    !> 行列構造体（実数型）
    type(monolis_mat_val_R) :: R
    !> 行列構造体（複素数型）
    type(monolis_mat_val_C) :: C
    !> 行列構造体（セパレート CSR 構造）
    type(monolis_mat_separated_CSR) :: SCSR
    !> 行列構造体（CSR 構造）
    type(monolis_mat_CSR) :: CSR
    !> 行列構造体（CSC 構造）
    type(monolis_mat_CSC) :: CSC
  end type monolis_mat

contains

  !> 行列構造体の初期化処理関数
  subroutine monolis_mat_initialize(monoMAT)
    implicit none
    !> 行列構造体
    type(monolis_mat) :: monoMAT

    monoMAT%N = 0
    monoMAT%NP = 0
    monoMAT%NDOF = 0

    call monolis_mat_initialize_val_R(monoMAT%R)
    call monolis_mat_initialize_val_C(monoMAT%C)
    call monolis_mat_initialize_SCSR(monoMAT%SCSR)
    call monolis_mat_initialize_CSR(monoMAT%CSR)
    call monolis_mat_initialize_CSC(monoMAT%CSC)
  end subroutine monolis_mat_initialize

  !> 行列構造体の初期化処理関数（実数型）
  subroutine monolis_mat_initialize_val_R(R)
    implicit none
    !> 行列構造体
    type(monolis_mat_val_R) :: R

    call monolis_dealloc_R_1d(R%A)
    call monolis_dealloc_R_1d(R%U)
    call monolis_dealloc_R_1d(R%D)
    call monolis_dealloc_R_1d(R%L)
    call monolis_dealloc_R_1d(R%X)
    call monolis_dealloc_R_1d(R%B)
  end subroutine monolis_mat_initialize_val_R

  !> 行列構造体の初期化処理関数（複素数型）
  subroutine monolis_mat_initialize_val_C(C)
    implicit none
    !> 行列構造体
    type(monolis_mat_val_C) :: C

    call monolis_dealloc_C_1d(C%A)
    call monolis_dealloc_C_1d(C%U)
    call monolis_dealloc_C_1d(C%D)
    call monolis_dealloc_C_1d(C%L)
    call monolis_dealloc_C_1d(C%X)
    call monolis_dealloc_C_1d(C%B)
  end subroutine monolis_mat_initialize_val_C

  !> 行列構造体の初期化処理関数（セパレート CSR 構造）
  subroutine monolis_mat_initialize_SCSR(SCSR)
    implicit none
    !> 行列構造体
    type(monolis_mat_separated_CSR) :: SCSR

    call monolis_dealloc_I_1d(SCSR%indexU)
    call monolis_dealloc_I_1d(SCSR%itemU)
    call monolis_dealloc_I_1d(SCSR%indexL)
    call monolis_dealloc_I_1d(SCSR%itemL)
  end subroutine monolis_mat_initialize_SCSR

  !> 行列構造体の初期化処理関数（CSR 構造）
  subroutine monolis_mat_initialize_CSR(CSR)
    implicit none
    !> 行列構造体
    type(monolis_mat_CSR) :: CSR

    call monolis_dealloc_I_1d(CSR%index)
    call monolis_dealloc_I_1d(CSR%item)
  end subroutine monolis_mat_initialize_CSR

  !> 行列構造体の初期化処理関数（CSC 構造）
  subroutine monolis_mat_initialize_CSC(CSC)
    implicit none
    !> 行列構造体
    type(monolis_mat_CSC) :: CSC

    call monolis_dealloc_I_1d(CSC%index)
    call monolis_dealloc_I_1d(CSC%item)
    call monolis_dealloc_I_1d(CSC%perm)
  end subroutine monolis_mat_initialize_CSC

  !> 行列構造体の終了処理関数
  subroutine monolis_mat_finalize(monoMAT)
    implicit none
    !> 行列構造体
    type(monolis_mat) :: monoMAT

    monoMAT%N = 0
    monoMAT%NP = 0
    monoMAT%NDOF = 0

    call monolis_mat_finalize_val_R(monoMAT%R)
    call monolis_mat_finalize_val_C(monoMAT%C)
    call monolis_mat_finalize_SCSR(monoMAT%SCSR)
    call monolis_mat_finalize_CSR(monoMAT%CSR)
    call monolis_mat_finalize_CSC(monoMAT%CSC)
  end subroutine monolis_mat_finalize

  !> 行列構造体の終了処理関数（実数型）
  subroutine monolis_mat_finalize_val_R(R)
    implicit none
    !> 行列構造体
    type(monolis_mat_val_R) :: R

    call monolis_dealloc_R_1d(R%A)
    call monolis_dealloc_R_1d(R%U)
    call monolis_dealloc_R_1d(R%D)
    call monolis_dealloc_R_1d(R%L)
    call monolis_dealloc_R_1d(R%X)
    call monolis_dealloc_R_1d(R%B)
  end subroutine monolis_mat_finalize_val_R

  !> 行列構造体の終了処理関数（複素数型）
  subroutine monolis_mat_finalize_val_C(C)
    implicit none
    !> 行列構造体
    type(monolis_mat_val_C) :: C

    call monolis_dealloc_C_1d(C%A)
    call monolis_dealloc_C_1d(C%U)
    call monolis_dealloc_C_1d(C%D)
    call monolis_dealloc_C_1d(C%L)
    call monolis_dealloc_C_1d(C%X)
    call monolis_dealloc_C_1d(C%B)
  end subroutine monolis_mat_finalize_val_C

  !> 行列構造体の終了処理関数（セパレート CSR 構造）
  subroutine monolis_mat_finalize_SCSR(SCSR)
    implicit none
    !> 行列構造体
    type(monolis_mat_separated_CSR) :: SCSR

    call monolis_dealloc_I_1d(SCSR%indexU)
    call monolis_dealloc_I_1d(SCSR%itemU)
    call monolis_dealloc_I_1d(SCSR%indexL)
    call monolis_dealloc_I_1d(SCSR%itemL)
  end subroutine monolis_mat_finalize_SCSR

  !> 行列構造体の終了処理関数（CSR 構造）
  subroutine monolis_mat_finalize_CSR(CSR)
    implicit none
    !> 行列構造体
    type(monolis_mat_CSR) :: CSR

    call monolis_dealloc_I_1d(CSR%index)
    call monolis_dealloc_I_1d(CSR%item)
  end subroutine monolis_mat_finalize_CSR

  !> 行列構造体の終了処理関数（CSC 構造）
  subroutine monolis_mat_finalize_CSC(CSC)
    implicit none
    !> 行列構造体
    type(monolis_mat_CSC) :: CSC

    call monolis_dealloc_I_1d(CSC%index)
    call monolis_dealloc_I_1d(CSC%item)
    call monolis_dealloc_I_1d(CSC%perm)
  end subroutine monolis_mat_finalize_CSC

  !> 右辺ベクトルの設定（実数型）
  subroutine monolis_set_RHS_R(monoMAT, B)
    implicit none
    !> 疎行列構造体
    type(monolis_mat) :: monoMAT
    !> 右辺ベクトル
    real(kdouble) :: B(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_set_RHS_R")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%R%B(i) = B(i)
    enddo
  end subroutine monolis_set_RHS_R

  !> 右辺ベクトルの設定（複素数型）
  subroutine monolis_set_RHS_C(monoMAT, B)
    implicit none
    !> 疎行列構造体
    type(monolis_mat) :: monoMAT
    !> 右辺ベクトル
    complex(kdouble) :: B(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_set_RHS_C")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%C%B(i) = B(i)
    enddo
  end subroutine monolis_set_RHS_C

  !> 解ベクトルの初期値設定（実数型）
  subroutine monolis_set_initial_solution_R(monoMAT, X)
    implicit none
    !> 疎行列構造体
    type(monolis_mat) :: monoMAT
    !> 解ベクトル
    real(kdouble) :: X(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_set_initial_solution_R")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%R%X(i) = X(i)
    enddo
  end subroutine monolis_set_initial_solution_R

  !> 解ベクトルの初期値設定（複素数型）
  subroutine monolis_set_initial_solution_C(monoMAT, X)
    implicit none
    !> 疎行列構造体
    type(monolis_mat) :: monoMAT
    !> 解ベクトル
    complex(kdouble) :: X(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_set_initial_solution_C")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%C%X(i) = X(i)
    enddo
  end subroutine monolis_set_initial_solution_C

  !> 解ベクトルの設定（実数型）
  subroutine monolis_get_solution_R(monoMAT, X)
    implicit none
    !> 疎行列構造体
    type(monolis_mat) :: monoMAT
    !> 解ベクトル
    real(kdouble) :: X(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_get_solution_R")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      X(i) = monoMAT%R%X(i)
    enddo
  end subroutine monolis_get_solution_R

  !> 解ベクトルの設定（複素数型）
  subroutine monolis_get_solution_C(monoMAT, X)
    implicit none
    !> 疎行列構造体
    type(monolis_mat) :: monoMAT
    !> 解ベクトル
    complex(kdouble) :: X(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_get_solution_C")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      X(i) = monoMAT%C%X(i)
    enddo
  end subroutine monolis_get_solution_C

  !> 入力パラメータのチェック
  subroutine monolis_check_input_param(monoCOM, monoMAT)
    implicit none
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT

    call monolis_std_debug_log_header("monolis_check_input_param")

    !if(monolis_mpi_local_comm_size(monoCOM%comm) > 1)then
    !endif

    if(monoMAT%N == 0)then
      call monolis_std_warning_string("the number of internal vertex is 0")
    endif
  end subroutine monolis_check_input_param

end module mod_monolis_def_mat
