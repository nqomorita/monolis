!> 行列構造体の定義
module mod_monolis_def_mat
  use mod_monolis_utils
  implicit none

#ifdef WITH_MUMPS
  include 'dmumps_struc.h'
#endif

  !> 行列構造体（実数型）
  type monolis_mat_val_R
    !> 全行列値
    real(kdouble), pointer, contiguous :: A(:) => null()
    !> 行列値（狭義上三角）
    real(kdouble), pointer, contiguous :: U(:) => null()
    !> 行列値（対角成分）
    real(kdouble), pointer, contiguous :: D(:) => null()
    !> 行列値（狭義下三角）
    real(kdouble), pointer, contiguous :: L(:) => null()
    !> 解ベクトル
    real(kdouble), pointer, contiguous :: X(:) => null()
    !> 右辺ベクトル
    real(kdouble), pointer, contiguous :: B(:) => null()
    !> 行列値（DIA 形式、column-major で (N, Ndiag) を 1 次元化）
    real(kdouble), pointer, contiguous :: Adia(:) => null()
    !> 行列値（ELL 形式、column-major で (N, Nmaxcol) を 1 次元化）
    real(kdouble), pointer, contiguous :: Aell(:) => null()
  end type monolis_mat_val_R

  !> 行列構造体（複素数型）
  type monolis_mat_val_C
    !> 全行列値
    complex(kdouble), pointer, contiguous :: A(:) => null()
    !> 行列値（狭義上三角）
    complex(kdouble), pointer, contiguous :: U(:) => null()
    !> 行列値（対角成分）
    complex(kdouble), pointer, contiguous :: D(:) => null()
    !> 行列値（狭義下三角）
    complex(kdouble), pointer, contiguous :: L(:) => null()
    !> 解ベクトル
    complex(kdouble), pointer, contiguous :: X(:) => null()
    !> 右辺ベクトル
    complex(kdouble), pointer, contiguous :: B(:) => null()
    !> 行列値（DIA 形式、column-major で (N, Ndiag) を 1 次元化）
    complex(kdouble), pointer, contiguous :: Adia(:) => null()
    !> 行列値（ELL 形式、column-major で (N, Nmaxcol) を 1 次元化）
    complex(kdouble), pointer, contiguous :: Aell(:) => null()
  end type monolis_mat_val_C

  !> 行列構造体（セパレート CSR 構造）
  type monolis_mat_separated_CSR
    !> index 配列（狭義上三角）
    integer(kint), pointer, contiguous :: indexU(:) => null()
    !> item 配列（狭義上三角）
    integer(kint), pointer, contiguous :: itemU(:) => null()
    !> index 配列（狭義下三角）
    integer(kint), pointer, contiguous :: indexL(:) => null()
    !> item 配列（狭義下三角）
    integer(kint), pointer, contiguous :: itemL(:) => null()
  end type monolis_mat_separated_CSR

  !> 行列構造体（CSR 構造）
  type monolis_mat_CSR
    !> index 配列
    integer(kint), pointer, contiguous :: index(:) => null()
    !> item 配列
    integer(kint), pointer, contiguous :: item(:) => null()
  end type monolis_mat_CSR

  !> 行列構造体（DIA 構造、整数情報）
  type monolis_mat_DIA
    !> 対角線本数
    integer(kint) :: Ndiag = 0
    !> 各対角の主対角からのオフセット（負=下三角、正=上三角）
    integer(kint), pointer, contiguous :: offset(:) => null()
  end type monolis_mat_DIA

  !> 行列構造体（ELL 構造、整数情報）
  type monolis_mat_ELL
    !> 1 行あたりの最大非ゼロブロック数
    integer(kint) :: Nmaxcol = 0
    !> 各スロットのブロック列番号（column-major で (N, Nmaxcol) を 1 次元化、0=パディング）
    integer(kint), pointer, contiguous :: col(:) => null()
  end type monolis_mat_ELL

  !> 行列構造体（CSC 構造）
  type monolis_mat_CSC
    !> index 配列
    integer(kint), pointer, contiguous :: index(:) => null()
    !> item 配列
    integer(kint), pointer, contiguous :: item(:) => null()
    !> CSR 形式に対する行列値の置換ベクトル
    integer(kint), pointer, contiguous :: perm(:) => null()
  end type monolis_mat_CSC

  !> 行列構造体（reordering 構造）
  type monolis_mat_reorder
    !> perm 配列
    integer(kint), pointer, contiguous :: perm(:) => null()
    !> iperm 配列
    integer(kint), pointer, contiguous :: iperm(:) => null()
    !> rperm 配列
    real(kdouble), pointer, contiguous :: rperm(:) => null()
    !> cperm 配列
    complex(kdouble), pointer, contiguous :: cperm(:) => null()
  end type monolis_mat_reorder

  type monolis_mat_DMUMPS
    integer(kint), allocatable :: offset_list(:)
    integer(kint), allocatable :: offset_counts(:)
    logical :: is_factored = .false.
    logical :: is_self = .false.
#ifdef WITH_MUMPS
    type (dmumps_struc), allocatable :: mumps(:)
#endif
  end type monolis_mat_DMUMPS

  !> 行列構造体（フロンタル行列構造）
  type :: monolis_mat_frontal
    !> フロント行列の全サイズ（ピボット + 更新）
    integer(kint) :: front_size = 0
    !> ピボット数
    integer(kint) :: pivot_size = 0
    !> 更新ブロックの行/列数
    integer(kint) :: update_size = 0
    !> ピボット列の LU 因子（front_size × pivot_size）
    real(kdouble), allocatable :: factor(:,:)
    !> ピボット行 × 更新列の U ブロック（pivot_size × update_size）
    real(kdouble), allocatable :: upper_update(:,:)
    !> 親フロントへの寄与（update_size × update_size、数値分解中の一時領域）
    real(kdouble), allocatable :: contribution(:,:)
  end type

  !> 行列構造体（LU 分解構造）
  type :: monolis_mat_lu
    !> 行列次元
    integer(kint) :: N = 0
    !> 置換配列（元 → 新、1-based）
    integer(kint), allocatable :: perm(:)
    !> 逆置換配列（新 → 元、1-based）
    integer(kint), allocatable :: iperm(:)
    !> フロント数
    integer(kint) :: nfronts = 0
    !> 最大フロントサイズ
    integer(kint) :: max_front_size = 0
    !> 各フロント先頭の置換後列番号
    integer(kint), allocatable :: super_start(:)
    !> フロント変数列 CSR の index（nfronts+1、1-based）
    integer(kint), allocatable :: front_ptr(:)
    !> フロント変数列（置換後の列番号、front_ptr で区切る）
    integer(kint), allocatable :: front_ind(:)
    !> フロント親
    integer(kint), allocatable :: front_parent(:)
    !> フロント第一子
    integer(kint), allocatable :: front_first_child(:)
    !> フロント次兄弟
    integer(kint), allocatable :: front_next_sibling(:)
    !> フロント後順序
    integer(kint), allocatable :: front_postorder(:)
    !> フロント全サイズ
    integer(kint), allocatable :: front_size(:)
    !> フロントごとのピボット数
    integer(kint), allocatable :: front_pivot_size(:)
    !> フロントごとの更新サイズ
    integer(kint), allocatable :: front_update_size(:)
    !> 元行列非ゼロ → フロント対応の index（nfronts+1、1-based）
    integer(kint), allocatable :: orig_ptr(:)
    !> 元行列非ゼロのフロント内行位置
    integer(kint), allocatable :: orig_row_pos(:)
    !> 元行列非ゼロのフロント内列位置
    integer(kint), allocatable :: orig_col_pos(:)
    !> 元行列非ゼロの A 配列内位置
    integer(kint), allocatable :: orig_entry(:)
    !> 子フロント寄与 → 親位置対応の index（nfronts+1、1-based）
    integer(kint), allocatable :: contrib_pos_ptr(:)
    !> 子の更新列が親で占める位置
    integer(kint), allocatable :: contrib_parent_pos(:)
    !> 子フロント連続区間 run の index（nfronts+1、1-based）
    integer(kint), allocatable :: contrib_run_ptr(:)
    !> run の子側開始位置
    integer(kint), allocatable :: contrib_run_first(:)
    !> run の長さ
    integer(kint), allocatable :: contrib_run_len(:)
    !> run の親側開始位置
    integer(kint), allocatable :: contrib_run_parent_first(:)
    !> 数値フロント配列
    type(monolis_mat_frontal), allocatable :: factors(:)
    !> 解析フェーズ完了フラグ
    logical :: analyzed = .false.
    !> 数値分解フェーズ完了フラグ
    logical :: factorized = .false.
  end type

  !> 行列構造体
  type monolis_mat
    !> 内部自由度数
    integer(kint) :: N
    !> 全自由度数
    integer(kint) :: NP
    !> 1 ブロックの自由度
    integer(kint) :: NDOF
    !> 1 ブロックの自由度配列
    integer(kint), pointer, contiguous :: n_dof_list(:) => null()
    !> 1 ブロックの自由度配列（index 型の圧縮形式）
    integer(kint), pointer, contiguous :: n_dof_index(:) => null()
    !> 1 ブロックの自由度配列（index 型の圧縮形式、ブロック自由度の 2 乗値）
    integer(kint), pointer, contiguous :: n_dof_index2(:) => null()
    !> 行列構造体（実数型）
    type(monolis_mat_val_R) :: R
    !> 行列構造体（複素数型）
    type(monolis_mat_val_C) :: C
    !> 行列構造体（セパレート CSR 構造）
    type(monolis_mat_separated_CSR) :: SCSR
    !> 行列構造体（CSR 構造）
    type(monolis_mat_CSR) :: CSR
    !> 行列構造体（DIA 構造、整数情報）
    type(monolis_mat_DIA) :: DIA
    !> 行列構造体（ELL 構造、整数情報）
    type(monolis_mat_ELL) :: ELL
    !> 行列構造体（CSC 構造）
    type(monolis_mat_CSC) :: CSC
    !> 行列構造体（reordering 構造）
    type(monolis_mat_reorder) :: REORDER
    !> 行列構造体（DMUMPS 構造）
    type(monolis_mat_DMUMPS) :: DMUMPS
    !> 行列構造体（LU 分解構造）
    type(monolis_mat_lu) :: LU
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

    call monolis_pdealloc_I_1d(monoMAT%n_dof_list)
    call monolis_pdealloc_I_1d(monoMAT%n_dof_index)
    call monolis_pdealloc_I_1d(monoMAT%n_dof_index2)

    call monolis_mat_initialize_val_R(monoMAT%R)
    call monolis_mat_initialize_val_C(monoMAT%C)
    call monolis_mat_initialize_SCSR(monoMAT%SCSR)
    call monolis_mat_initialize_CSR(monoMAT%CSR)
    call monolis_mat_initialize_DIA(monoMAT%DIA)
    call monolis_mat_initialize_ELL(monoMAT%ELL)
    call monolis_mat_initialize_CSC(monoMAT%CSC)
    call monolis_mat_initialize_REORDER(monoMAT%REORDER)
    call monolis_mat_initialize_LU(monoMAT%LU)
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
    call monolis_pdealloc_R_1d(R%Adia)
    call monolis_pdealloc_R_1d(R%Aell)
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
    call monolis_pdealloc_C_1d(C%Adia)
    call monolis_pdealloc_C_1d(C%Aell)
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
  !> 行列構造体の初期化処理関数（DIA 構造）
  subroutine monolis_mat_initialize_DIA(DIA)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_DIA), intent(inout) :: DIA

    DIA%Ndiag = 0
    call monolis_pdealloc_I_1d(DIA%offset)
  end subroutine monolis_mat_initialize_DIA

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数（DIA 構造）
  subroutine monolis_mat_finalize_DIA(DIA)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_DIA), intent(inout) :: DIA

    DIA%Ndiag = 0
    call monolis_pdealloc_I_1d(DIA%offset)
  end subroutine monolis_mat_finalize_DIA

  !> @ingroup def_mat_init
  !> 行列構造体の初期化処理関数（ELL 構造）
  subroutine monolis_mat_initialize_ELL(ELL)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_ELL), intent(inout) :: ELL

    ELL%Nmaxcol = 0
    call monolis_pdealloc_I_1d(ELL%col)
  end subroutine monolis_mat_initialize_ELL

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数（ELL 構造）
  subroutine monolis_mat_finalize_ELL(ELL)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_ELL), intent(inout) :: ELL

    ELL%Nmaxcol = 0
    call monolis_pdealloc_I_1d(ELL%col)
  end subroutine monolis_mat_finalize_ELL

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
  !> 行列構造体の初期化処理関数（REORDER 構造）
  subroutine monolis_mat_initialize_REORDER(REORDER)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_REORDER), intent(inout) :: REORDER

    call monolis_pdealloc_I_1d(REORDER%perm)
    call monolis_pdealloc_I_1d(REORDER%iperm)
    call monolis_pdealloc_R_1d(REORDER%rperm)
    call monolis_pdealloc_C_1d(REORDER%cperm)
  end subroutine monolis_mat_initialize_REORDER

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数
  subroutine monolis_mat_finalize(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT

    monoMAT%N = 0
    monoMAT%NP = 0
    monoMAT%NDOF = 0

    call monolis_pdealloc_I_1d(monoMAT%n_dof_list)
    call monolis_pdealloc_I_1d(monoMAT%n_dof_index)
    call monolis_pdealloc_I_1d(monoMAT%n_dof_index2)

    call monolis_mat_finalize_val_R(monoMAT%R)
    call monolis_mat_finalize_val_C(monoMAT%C)
    call monolis_mat_finalize_SCSR(monoMAT%SCSR)
    call monolis_mat_finalize_CSR(monoMAT%CSR)
    call monolis_mat_finalize_DIA(monoMAT%DIA)
    call monolis_mat_finalize_ELL(monoMAT%ELL)
    call monolis_mat_finalize_CSC(monoMAT%CSC)
    call monolis_mat_finalize_REORDER(monoMAT%REORDER)
    call monolis_mat_finalize_MUMPS(monoMAT%DMUMPS)
    call monolis_mat_finalize_LU(monoMAT%LU)
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
    call monolis_pdealloc_R_1d(R%Adia)
    call monolis_pdealloc_R_1d(R%Aell)
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
    call monolis_pdealloc_C_1d(C%Adia)
    call monolis_pdealloc_C_1d(C%Aell)
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
  !> 行列構造体の初期化処理関数（REORDER 構造）
  subroutine monolis_mat_finalize_REORDER(REORDER)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_REORDER), intent(inout) :: REORDER

    call monolis_pdealloc_I_1d(REORDER%perm)
    call monolis_pdealloc_I_1d(REORDER%iperm)
    call monolis_pdealloc_R_1d(REORDER%rperm)
    call monolis_pdealloc_C_1d(REORDER%cperm)
  end subroutine monolis_mat_finalize_REORDER

  !> @ingroup def_mat_init
  !> 行列構造体の初期化処理関数（MUMPS 構造）
  subroutine monolis_mat_finalize_MUMPS(MUMPS)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_DMUMPS), intent(inout) :: MUMPS

    call monolis_dealloc_I_1d(MUMPS%offset_list)
    call monolis_dealloc_I_1d(MUMPS%offset_counts)

#ifdef WITH_MUMPS
    if(allocated(MUMPS%mumps))then
      MUMPS%mumps(1)%JOB = -2
      call DMUMPS(MUMPS%mumps(1))
    endif
#endif
  end subroutine monolis_mat_finalize_MUMPS

  !> @ingroup def_mat_init
  !> 行列構造体の初期化処理関数（フロンタル行列構造）
  subroutine monolis_mat_initialize_frontal(FRONTAL)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_frontal), intent(inout) :: FRONTAL

    FRONTAL%front_size = 0
    FRONTAL%pivot_size = 0
    FRONTAL%update_size = 0
    call monolis_dealloc_R_2d(FRONTAL%factor)
    call monolis_dealloc_R_2d(FRONTAL%upper_update)
    call monolis_dealloc_R_2d(FRONTAL%contribution)
  end subroutine monolis_mat_initialize_frontal

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数（フロンタル行列構造）
  subroutine monolis_mat_finalize_frontal(FRONTAL)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_frontal), intent(inout) :: FRONTAL

    FRONTAL%front_size = 0
    FRONTAL%pivot_size = 0
    FRONTAL%update_size = 0
    call monolis_dealloc_R_2d(FRONTAL%factor)
    call monolis_dealloc_R_2d(FRONTAL%upper_update)
    call monolis_dealloc_R_2d(FRONTAL%contribution)
  end subroutine monolis_mat_finalize_frontal

  !> @ingroup def_mat_init
  !> 行列構造体の初期化処理関数（LU 分解構造）
  subroutine monolis_mat_initialize_LU(LU)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_lu), intent(inout) :: LU

    LU%N = 0
    LU%nfronts = 0
    LU%max_front_size = 0
    LU%analyzed = .false.
    LU%factorized = .false.

    call monolis_mat_finalize_LU(LU)
  end subroutine monolis_mat_initialize_LU

  !> @ingroup def_mat_init
  !> 行列構造体の終了処理関数（LU 分解構造）
  subroutine monolis_mat_finalize_LU(LU)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat_lu), intent(inout) :: LU
    integer(kint) :: i

    if(allocated(LU%factors))then
      do i = 1, size(LU%factors)
        call monolis_mat_finalize_frontal(LU%factors(i))
      enddo
      deallocate(LU%factors)
    endif

    call monolis_dealloc_I_1d(LU%perm)
    call monolis_dealloc_I_1d(LU%iperm)
    call monolis_dealloc_I_1d(LU%super_start)
    call monolis_dealloc_I_1d(LU%front_ptr)
    call monolis_dealloc_I_1d(LU%front_ind)
    call monolis_dealloc_I_1d(LU%front_parent)
    call monolis_dealloc_I_1d(LU%front_first_child)
    call monolis_dealloc_I_1d(LU%front_next_sibling)
    call monolis_dealloc_I_1d(LU%front_postorder)
    call monolis_dealloc_I_1d(LU%front_size)
    call monolis_dealloc_I_1d(LU%front_pivot_size)
    call monolis_dealloc_I_1d(LU%front_update_size)
    call monolis_dealloc_I_1d(LU%orig_ptr)
    call monolis_dealloc_I_1d(LU%orig_row_pos)
    call monolis_dealloc_I_1d(LU%orig_col_pos)
    call monolis_dealloc_I_1d(LU%orig_entry)
    call monolis_dealloc_I_1d(LU%contrib_pos_ptr)
    call monolis_dealloc_I_1d(LU%contrib_parent_pos)
    call monolis_dealloc_I_1d(LU%contrib_run_ptr)
    call monolis_dealloc_I_1d(LU%contrib_run_first)
    call monolis_dealloc_I_1d(LU%contrib_run_len)
    call monolis_dealloc_I_1d(LU%contrib_run_parent_first)

    LU%analyzed = .false.
    LU%factorized = .false.
    LU%nfronts = 0
    LU%max_front_size = 0
    LU%N = 0
  end subroutine monolis_mat_finalize_LU

  !> @ingroup def_mat_init
  !> 右辺ベクトルの設定（実数型）
  subroutine monolis_set_RHS_R(monoMAT, B)
    implicit none
    !> [in,out] 疎行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: B(:)
    integer(kint) :: i, NNDOF, NPNDOF

    call monolis_std_debug_log_header("monolis_set_RHS_R")

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    do i = 1, NPNDOF
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
    integer(kint) :: i, NNDOF, NPNDOF

    call monolis_std_debug_log_header("monolis_set_RHS_C")

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    do i = 1, NPNDOF
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
    integer(kint) :: i, NNDOF, NPNDOF

    call monolis_std_debug_log_header("monolis_set_initial_solution_R")

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    do i = 1, NPNDOF
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
    integer(kint) :: i, NNDOF, NPNDOF

    call monolis_std_debug_log_header("monolis_set_initial_solution_C")

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    do i = 1, NPNDOF
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
    integer(kint) :: i, NNDOF, NPNDOF

    call monolis_std_debug_log_header("monolis_get_solution_R")

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    do i = 1, NPNDOF
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
    integer(kint) :: i, NNDOF, NPNDOF

    call monolis_std_debug_log_header("monolis_get_solution_C")

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    do i = 1, NPNDOF
      X(i) = monoMAT%C%X(i)
    enddo
  end subroutine monolis_get_solution_C

  !> @ingroup def_mat_init
  !> 疎行列に対応するベクトル配列サイズの取得関数
  subroutine monolis_get_vec_size(N, NP, NDOF, n_dof_index, N_size, NP_size)
    implicit none
    !> [in] 内部計算点数
    integer(kint), intent(in) :: N
    !> [in] 全計算点数
    integer(kint), intent(in) :: NP
    !> [in] 自由度（固定値）
    integer(kint), intent(in) :: NDOF
    !> [in] 自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [out] ベクトルサイズ（内部計算点数相当）
    integer(kint), intent(out) :: N_size
    !> [out] ベクトルサイズ（全計算点数相当）
    integer(kint), intent(out) :: NP_size

    if(NDOF == -1)then
      n_size  = n_dof_index(N  + 1)
      np_size = n_dof_index(NP + 1)
    else
      n_size  = N *NDOF
      np_size = NP*NDOF
    endif
  end subroutine monolis_get_vec_size
  !> @ingroup def_mat_init

  !> 疎行列に対応するベクトル配列サイズの取得関数
  subroutine monolis_get_mat_size(MAT, NZ, NZD, NZU, NZL)
    implicit none
    !> [in,out] 疎行列構造体
    type(monolis_mat), intent(in) :: MAT
    !> [in] 内部計算点数
    integer(kint), intent(out) :: NZ
    !> [in] 全計算点数
    integer(kint), intent(out) :: NZD
    !> [in] 自由度（固定値）
    integer(kint), intent(out) :: NZU
    !> [in] 自由度（固定値）
    integer(kint), intent(out) :: NZL
    integer(kint) :: i, j, in, jS, jE, NP, NDOF

    NP = MAT%NP
    NDOF = MAT%NDOF

    if(NDOF == -1)then
      in = MAT%CSR%index(NP + 1)
      NZ = MAT%n_dof_index2(in + 1)

      NZD = 0
      NZU = 0
      NZL = 0
      do i = 1, NP
        jS = MAT%CSR%index(i) + 1
        jE = MAT%CSR%index(i + 1)
        do j = jS, jE
          in = MAT%CSR%item(j)
          if(i == in)then
            NZD = NZD + MAT%n_dof_index(i)*MAT%n_dof_index(in)
          elseif(i < in)then
            NZU = NZU + MAT%n_dof_index(i)*MAT%n_dof_index(in)
          elseif(in < i)then
            NZL = NZL + MAT%n_dof_index(i)*MAT%n_dof_index(in)
          endif
        enddo
      enddo

    else
      NZ = MAT%CSR%index(NP + 1)
      NZ = NZ*NDOF*NDOF

      NZU = 0
      if(associated(MAT%SCSR%indexU))then
        NZU = MAT%SCSR%indexU(NP + 1)
      endif
  
      NZL = 0
      if(associated(MAT%SCSR%indexL))then
        NZL = MAT%SCSR%indexL(NP + 1)
      endif

      NZD = NP*NDOF*NDOF
      NZU = NZU*NDOF*NDOF
      NZL = NZL*NDOF*NDOF
    endif
  end subroutine monolis_get_mat_size
end module mod_monolis_def_mat
