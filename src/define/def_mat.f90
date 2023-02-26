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
    !> 内部自由度数
    integer(kint) :: N
    !> 全自由度数
    integer(kint) :: NP
    !> 非零個数（狭義上三角）
    integer(kint) :: NPU
    !> 非零個数（狭義下三角）
    integer(kint) :: NPL
    !> 1 ブロックの自由度
    integer(kint) :: NDOF
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
    !> 内部自由度数
    integer(kint) :: N
    !> 全自由度数
    integer(kint) :: NP
    !> 非零個数
    integer(kint) :: NPZ
    !> 1 ブロックの自由度
    integer(kint) :: NDOF
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
    type(monolis_mat_val_R) :: R
    type(monolis_mat_val_C) :: C
    type(monolis_mat_separated_CSR) :: SCSR
    type(monolis_mat_CSR) :: CSR
    type(monolis_mat_CSC) :: CSC
  end type monolis_mat

  !> 前処理行列構造体
  type monolis_prec
    type(monolis_mat_val_R) :: R
    type(monolis_mat_val_C) :: C
    type(monolis_mat_separated_CSR) :: SCSR
    type(monolis_mat_CSR) :: CSR
    type(monolis_mat_CSC) :: CSC
  end type monolis_prec

contains

  !> 行列構造体の初期化処理関数
  subroutine monolis_mat_initialize(monoMAT)
    implicit none
    !> 行列構造体
    type(monolis_mat) :: monoMAT

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

    SCSR%N = 0
    SCSR%NP = 0
    SCSR%NPU = 0
    SCSR%NPL = 0
    SCSR%NDOF = 0
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

    CSR%N = 0
    CSR%NP = 0
    CSR%NPZ = 0
    CSR%NDOF = 0
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

    SCSR%N = 0
    SCSR%NP = 0
    SCSR%NPU = 0
    SCSR%NPL = 0
    SCSR%NDOF = 0
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

    CSR%N = 0
    CSR%NP = 0
    CSR%NPZ = 0
    CSR%NDOF = 0
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
end module mod_monolis_def_mat