!> 行列構造体の定義
module mod_monolis_def_mat
  use mod_monolis_utils
  implicit none

  !> 行列構造体（実数型）
  type monolis_mat_val_R
    !> 全行列値
    real(kdouble), pointer :: A(:) => null()
    !> 行列値（狭義上三角）
    real(kdouble), pointer :: U(:) => null()
    !> 行列値（対角成分）
    real(kdouble), pointer :: D(:) => null()
    !> 行列値（狭義下三角）
    real(kdouble), pointer :: L(:) => null()
    !> 解ベクトル
    real(kdouble), pointer :: X(:) => null()
    !> 右辺ベクトル
    real(kdouble), pointer :: B(:) => null()
  end type monolis_mat_val_R

  !> 行列構造体（複素数型）
  type monolis_mat_val_C
    !> 全行列値
    complex(kdouble), pointer :: A(:) => null()
    !> 行列値（狭義上三角）
    complex(kdouble), pointer :: U(:) => null()
    !> 行列値（対角成分）
    complex(kdouble), pointer :: D(:) => null()
    !> 行列値（狭義下三角）
    complex(kdouble), pointer :: L(:) => null()
    !> 解ベクトル
    complex(kdouble), pointer :: X(:) => null()
    !> 右辺ベクトル
    complex(kdouble), pointer :: B(:) => null()
  end type monolis_mat_val_C

  !> 行列構造体（セパレート CSR 構造）
  type monolis_mat_separated_CSR
    !> index 配列（狭義上三角）
    integer(kint), pointer :: indexU(:) => null()
    !> item 配列（狭義上三角）
    integer(kint), pointer :: itemU(:) => null()
    !> index 配列（狭義下三角）
    integer(kint), pointer :: indexL(:) => null()
    !> item 配列（狭義下三角）
    integer(kint), pointer :: itemL(:) => null()
  end type monolis_mat_separated_CSR

  !> 行列構造体（CSR 構造）
  type monolis_mat_CSR
    !> index 配列
    integer(kint), pointer :: index(:) => null()
    !> item 配列
    integer(kint), pointer :: item(:) => null()
  end type monolis_mat_CSR

  !> 行列構造体（CSC 構造）
  type monolis_mat_CSC
    !> index 配列
    integer(kint), pointer :: index(:) => null()
    !> item 配列
    integer(kint), pointer :: item(:) => null()
    !> CSR 形式に対する行列値の置換ベクトル
    integer(kint), pointer :: perm(:) => null()
  end type monolis_mat_CSC

  !> 行列構造体（格子構造）
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

  !> @ingroup def_init
  !> 行列構造体の初期化処理関数
  subroutine monolis_mat_initialize(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT

    monoMAT%N = 0
    monoMAT%NP = 0
    monoMAT%NDOF = 0

    call monolis_mat_initialize_val_R(monoMAT%R)
    call monolis_mat_initialize_val_C(monoMAT%C)
    call monolis_mat_initialize_SCSR(monoMAT%SCSR)
    call monolis_mat_initialize_CSR(monoMAT%CSR)
    call monolis_mat_initialize_CSC(monoMAT%CSC)
  end subroutine monolis_mat_initialize

  !> @ingroup def_mat_init
  !> 行列構造体の初期化処理関数（実数型）
  subroutine monolis_mat_initialize_val_R(R)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_val_R), intent(inout) :: R

    call monolis_pdealloc_R_1d(R%A)
    call monolis_pdealloc_R_1d(R%U)
    call monolis_pdealloc_R_1d(R%D)
    call monolis_pdealloc_R_1d(R%L)
    call monolis_pdealloc_R_1d(R%X)
    call monolis_pdealloc_R_1d(R%B)
  end subroutine monolis_mat_initialize_val_R

  !> @ingroup def_mat_init
  !> 行列構造体の初期化処理関数（複素数型）
  subroutine monolis_mat_initialize_val_C(C)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_val_C), intent(inout) :: C

    call monolis_pdealloc_C_1d(C%A)
    call monolis_pdealloc_C_1d(C%U)
    call monolis_pdealloc_C_1d(C%D)
    call monolis_pdealloc_C_1d(C%L)
    call monolis_pdealloc_C_1d(C%X)
    call monolis_pdealloc_C_1d(C%B)
  end subroutine monolis_mat_initialize_val_C

  !> @ingroup def_mat_init
  !> 行列構造体の初期化処理関数（セパレート CSR 構造）
  subroutine monolis_mat_initialize_SCSR(SCSR)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_separated_CSR), intent(inout) :: SCSR

    call monolis_pdealloc_I_1d(SCSR%indexU)
    call monolis_pdealloc_I_1d(SCSR%itemU)
    call monolis_pdealloc_I_1d(SCSR%indexL)
    call monolis_pdealloc_I_1d(SCSR%itemL)
  end subroutine monolis_mat_initialize_SCSR

  !> @ingroup def_mat_init
  !> 行列構造体の初期化処理関数（CSR 構造）
  subroutine monolis_mat_initialize_CSR(CSR)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_CSR), intent(inout) :: CSR

    call monolis_pdealloc_I_1d(CSR%index)
    call monolis_pdealloc_I_1d(CSR%item)
  end subroutine monolis_mat_initialize_CSR

  !> @ingroup def_mat_init
  !> 行列構造体の初期化処理関数（CSC 構造）
  subroutine monolis_mat_initialize_CSC(CSC)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_CSC), intent(inout) :: CSC

    call monolis_pdealloc_I_1d(CSC%index)
    call monolis_pdealloc_I_1d(CSC%item)
    call monolis_pdealloc_I_1d(CSC%perm)
  end subroutine monolis_mat_initialize_CSC

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数
  subroutine monolis_mat_finalize(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT

    monoMAT%N = 0
    monoMAT%NP = 0
    monoMAT%NDOF = 0

    call monolis_mat_finalize_val_R(monoMAT%R)
    call monolis_mat_finalize_val_C(monoMAT%C)
    call monolis_mat_finalize_SCSR(monoMAT%SCSR)
    call monolis_mat_finalize_CSR(monoMAT%CSR)
    call monolis_mat_finalize_CSC(monoMAT%CSC)
  end subroutine monolis_mat_finalize

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数（実数型）
  subroutine monolis_mat_finalize_val_R(R)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_val_R), intent(inout) :: R

    call monolis_pdealloc_R_1d(R%A)
    call monolis_pdealloc_R_1d(R%U)
    call monolis_pdealloc_R_1d(R%D)
    call monolis_pdealloc_R_1d(R%L)
    call monolis_pdealloc_R_1d(R%X)
    call monolis_pdealloc_R_1d(R%B)
  end subroutine monolis_mat_finalize_val_R

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数（複素数型）
  subroutine monolis_mat_finalize_val_C(C)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_val_C), intent(inout) :: C

    call monolis_pdealloc_C_1d(C%A)
    call monolis_pdealloc_C_1d(C%U)
    call monolis_pdealloc_C_1d(C%D)
    call monolis_pdealloc_C_1d(C%L)
    call monolis_pdealloc_C_1d(C%X)
    call monolis_pdealloc_C_1d(C%B)
  end subroutine monolis_mat_finalize_val_C

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数（セパレート CSR 構造）
  subroutine monolis_mat_finalize_SCSR(SCSR)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_separated_CSR), intent(inout) :: SCSR

    call monolis_pdealloc_I_1d(SCSR%indexU)
    call monolis_pdealloc_I_1d(SCSR%itemU)
    call monolis_pdealloc_I_1d(SCSR%indexL)
    call monolis_pdealloc_I_1d(SCSR%itemL)
  end subroutine monolis_mat_finalize_SCSR

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数（CSR 構造）
  subroutine monolis_mat_finalize_CSR(CSR)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_CSR), intent(inout) :: CSR

    call monolis_pdealloc_I_1d(CSR%index)
    call monolis_pdealloc_I_1d(CSR%item)
  end subroutine monolis_mat_finalize_CSR

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数（CSC 構造）
  subroutine monolis_mat_finalize_CSC(CSC)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_CSC), intent(inout) :: CSC

    call monolis_pdealloc_I_1d(CSC%index)
    call monolis_pdealloc_I_1d(CSC%item)
    call monolis_pdealloc_I_1d(CSC%perm)
  end subroutine monolis_mat_finalize_CSC

  !> @ingroup def_mat_init
  !> 右辺ベクトルの設定（実数型）
  subroutine monolis_set_RHS_R(monoMAT, B)
    implicit none
    !> [in,out] 疎行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: B(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_set_RHS_R")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%R%B(i) = B(i)
    enddo
  end subroutine monolis_set_RHS_R

  !> @ingroup def_mat_init
  !> 右辺ベクトルの設定（複素数型）
  subroutine monolis_set_RHS_C(monoMAT, B)
    implicit none
    !> [in,out] 疎行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in] 右辺ベクトル
    complex(kdouble), intent(in) :: B(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_set_RHS_C")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%C%B(i) = B(i)
    enddo
  end subroutine monolis_set_RHS_C

  !> @ingroup def_mat_init
  !> 解ベクトルの初期値設定（実数型）
  subroutine monolis_set_initial_solution_R(monoMAT, X)
    implicit none
    !> [in,out] 疎行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in] 解ベクトル
    real(kdouble), intent(in) :: X(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_set_initial_solution_R")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%R%X(i) = X(i)
    enddo
  end subroutine monolis_set_initial_solution_R

  !> @ingroup def_mat_init
  !> 解ベクトルの初期値設定（複素数型）
  subroutine monolis_set_initial_solution_C(monoMAT, X)
    implicit none
    !> [in,out] 疎行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in] 解ベクトル
    complex(kdouble), intent(in) :: X(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_set_initial_solution_C")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%C%X(i) = X(i)
    enddo
  end subroutine monolis_set_initial_solution_C

  !> @ingroup def_mat_init
  !> 解ベクトルの設定（実数型）
  subroutine monolis_get_solution_R(monoMAT, X)
    implicit none
    !> [in,out] 疎行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [out] 解ベクトル
    real(kdouble), intent(out) :: X(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_get_solution_R")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      X(i) = monoMAT%R%X(i)
    enddo
  end subroutine monolis_get_solution_R

  !> @ingroup def_mat_init
  !> 解ベクトルの設定（複素数型）
  subroutine monolis_get_solution_C(monoMAT, X)
    implicit none
    !> [in,out] 疎行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [out] 解ベクトル
    complex(kdouble), intent(out) :: X(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_get_solution_C")

    do i = 1, monoMAT%NP*monoMAT%NDOF
      X(i) = monoMAT%C%X(i)
    enddo
  end subroutine monolis_get_solution_C
end module mod_monolis_def_mat
