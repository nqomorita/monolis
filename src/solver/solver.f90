!> 線形ソルバモジュール
module mod_monolis_solve
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_solver_CG
  use mod_monolis_solver_BiCGSTAB
  use mod_monolis_solver_BiCGSTAB_N128
  use mod_monolis_solver_BiCGSTAB_noprec
  use mod_monolis_solver_GropCG
  use mod_monolis_solver_PipeCG
  use mod_monolis_solver_PipeCR
  use mod_monolis_solver_PipeBiCGSTAB
  use mod_monolis_solver_PipeBiCGSTAB_noprec
  use mod_monolis_solver_SOR
  use mod_monolis_solver_BiCGSAFE
  use mod_monolis_solver_IDRS
  use mod_monolis_solver_COCG
  use mod_monolis_precond
  use mod_monolis_spmat_convert_dia
  use mod_monolis_spmat_convert_ell
  use mod_monolis_spmat_handler_util
#ifdef _OPENACC
  use openacc
#endif

  implicit none

contains

  !> @ingroup solver
  !> 線形ソルバ関数（実数型）
  subroutine monolis_solve_R(monolis, monoCOM, B, X)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: B(:)
    !> [in,out] 解ベクトル
    real(kdouble), intent(inout) :: X(:)

    call monolis_set_RHS_R(monolis%MAT, B)

    call monolis_set_initial_solution_R(monolis%MAT, X)

    if(monoCOM%comm_size > 1) monolis%MAT%N = monoCOM%n_internal_vertex

    call monolis_solve_main_R(monolis%PRM, monoCOM, monolis%MAT, monolis%PREC)

    call monolis_get_solution_R(monolis%MAT, X)
  end subroutine monolis_solve_R

  !> @ingroup solver
  !> 線形ソルバ関数（複素数型）
  subroutine monolis_solve_C(monolis, monoCOM, B, X)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 右辺ベクトル
    complex(kdouble), intent(in) :: B(:)
    !> [in,out] 解ベクトル
    complex(kdouble), intent(inout) :: X(:)

    call monolis_set_RHS_C(monolis%MAT, B)

    call monolis_set_initial_solution_C(monolis%MAT, X)

    if(monoCOM%comm_size > 1) monolis%MAT%N = monoCOM%n_internal_vertex

    call monolis_solve_main_C(monolis%PRM, monoCOM, monolis%MAT, monolis%PREC)

    call monolis_get_solution_C(monolis%MAT, X)
  end subroutine monolis_solve_C

  !> @ingroup dev_solver
  !> 線形ソルバ関数（メイン関数）
  subroutine monolis_solve_main_R(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC
#ifdef _OPENACC
    integer(kint) :: NNDOF, NPNDOF, method
    real(kdouble), pointer, contiguous :: X(:), B(:), precD(:)
    integer(kint), pointer, contiguous :: matNdofList(:), matNdofIndex(:)
    logical :: is_ell, is_var, use_device, was_present
#endif

    call monolis_std_debug_log_header("monolis_solve_main_R")

#ifdef _OPENACC
    !# GPU 正本モード:組み立て中の非同期カーネルを待ち、デバイス組立の遅延エラーを検査する
    !# (前処理生成・形式変換は組み立て完了後の行列値を参照する)
    call monolis_check_matrix_assemble_error()
#endif

    call monolis_timer_initialize(monoPRM)

    !call monolis_check_input_param(monoCOM, monoMAT)

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT, monoPREC)

#ifdef _OPENACC
    !# GPU 実装のない解法（BiCGSTAB_N128 / SOR / IDRS）はデバイス常駐を行わずホスト（CSR 形式）で実行する
    method = monoPRM%Iarray(monolis_prm_I_method)
    use_device = .not. (method == monolis_iter_BiCGSTAB_N128 &
      & .or. method == monolis_iter_SOR &
      & .or. method == monolis_iter_IDRS)

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)
    X => monoMAT%R%X
    B => monoMAT%R%B
    was_present = acc_is_present(X)

    if(use_device)then
    !# 初期解のゼロクリア(デバイス常駐の場合はデバイス側もクリアする。ライトスルー)
    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
      if(was_present)then
!$acc kernels present(X)
        X = 0.0d0
!$acc end kernels
      endif
    endif

    !# CSR 形式から DIA / ELL 形式へ変換し、行列・解・右辺・前処理対角をデバイスに常駐させる
    !# (GPU 正本モードでは値の詰め替えはデバイス上で実行される)
    is_ell = (monoPRM%Iarray(monolis_prm_I_spmv_format) /= monolis_spmv_DIA)
    is_var = (monoMAT%NDOF == -1)
    if(is_ell)then
      if(is_var)then
        call monolis_convert_CSR_to_ELL_V_R(monoMAT)
      else
        call monolis_convert_CSR_to_ELL_R(monoMAT)
      endif
    else
      if(is_var)then
        call monolis_convert_CSR_to_DIA_V_R(monoMAT)
      else
        call monolis_convert_CSR_to_DIA_R(monoMAT)
      endif
    endif

    precD => monoPREC%R%D

    !# デバイス常駐(GPU 正本)の場合はデバイス側の X・B が最新のため転送しない。
    !# 非常駐(レガシー)経路ではホストの右辺・初期解をデバイスへ転送する
    !$acc enter data create(X(1:NPNDOF), B(1:NPNDOF))
    if(.not. was_present)then
      !$acc update device(X(1:NPNDOF), B(1:NPNDOF))
    endif
    !$acc enter data copyin(precD)
    if(is_ell)then
      !$acc enter data copyin(monoMAT%R%Aell, monoMAT%ELL%col)
    else
      !$acc enter data copyin(monoMAT%R%Adia, monoMAT%DIA%offset)
    endif

    !# 可変ブロック版は値オフセットと自由度配列もデバイスに常駐させる
    if(is_var)then
      matNdofList  => monoMAT%n_dof_list
      matNdofIndex => monoMAT%n_dof_index
      !$acc enter data copyin(matNdofList, matNdofIndex)
      if(is_ell)then
        !$acc enter data copyin(monoMAT%ELL%Vptr)
      else
        !$acc enter data copyin(monoMAT%DIA%Vptr)
      endif
    endif

    else
    !# GPU 未対応解法のホスト実行:X・B の常駐を一時解除して
    !# 全演算をホスト(非常駐)経路に揃える(A のホスト鏡像はライトスルーにより最新)
    if(was_present)then
      !$acc update self(X(1:NPNDOF), B(1:NPNDOF))
      !$acc exit data delete(X, B)
    endif
    endif
#endif

    call monolis_solver_select_R(monoPRM, monoCOM, monoMAT, monoPREC)

#ifdef _OPENACC
    if(use_device)then
    !# ソルバの事後条件：デバイス上の X は共有節点（袖）を含めて最新である
    !# ホスト同期が有効な場合のみ X をホストへ書き戻す(非常駐経路では常に書き戻す)
    if(monoPRM%Iarray(monolis_prm_I_is_host_sync) == monolis_I_true .or. .not. was_present)then
      !$acc update self(X(1:NPNDOF))
    endif
    if(is_var)then
      if(is_ell)then
        !$acc exit data delete(monoMAT%ELL%Vptr)
      else
        !$acc exit data delete(monoMAT%DIA%Vptr)
      endif
      !$acc exit data delete(matNdofList, matNdofIndex)
    endif
    !$acc exit data delete(precD)
    !$acc exit data delete(X(1:NPNDOF), B(1:NPNDOF))
    if(is_ell)then
      !$acc exit data delete(monoMAT%R%Aell, monoMAT%ELL%col)
      call monolis_dealloc_ELL_R(monoMAT)
    else
      !$acc exit data delete(monoMAT%R%Adia, monoMAT%DIA%offset)
      call monolis_dealloc_DIA_R(monoMAT)
    endif

    else
    !# ホスト実行解法の結果(X・B)をデバイスへ書き戻して常駐を復元する
    if(was_present)then
      !$acc enter data copyin(X, B)
    endif
    endif
#endif

    call monolis_precond_clear(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_timer_finalize(monoPRM, monoCOM)
  end subroutine monolis_solve_main_R

  !> @ingroup dev_solver
  !> 線形ソルバ関数（メイン関数）
  subroutine monolis_solve_main_C(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_solve_main_C")

    call monolis_timer_initialize(monoPRM)

    !call monolis_check_input_param(monoCOM, monoMAT)

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_solver_select_C(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_precond_clear(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_timer_finalize(monoPRM, monoCOM)
  end subroutine monolis_solve_main_C

  !> @ingroup dev_solver
  !> 線形ソルバの選択関数
  subroutine monolis_solver_select_R(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_solver_select_R")

    if(monoPRM%Iarray(monolis_prm_I_show_summary) == monolis_I_true .and. monoCOM%my_rank == 0)then
      write(*,"(a)") &
      & "** MONOLIS solver: "//trim(monolis_str_iter(monoPRM%Iarray(monolis_prm_I_method)))//&
      & ", prec: "//trim(monolis_str_prec(monoPRM%Iarray(monolis_prm_I_precond)))
    endif

    select case(monoPRM%Iarray(monolis_prm_I_method))
      case (monolis_iter_CG)
        call monolis_solver_CG(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_BiCGSTAB)
#ifdef _OPENACC
        call monolis_solver_BiCGSTAB_GPU(monoPRM, monoCOM, monoMAT, monoPREC)
#else
        call monolis_solver_BiCGSTAB(monoPRM, monoCOM, monoMAT, monoPREC)
#endif

      case (monolis_iter_BiCGSTAB_noprec)
        call monolis_solver_BiCGSTAB_noprec(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_GropCG)
        call monolis_solver_GropCG(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_PipeCG)
        call monolis_solver_PipeCG(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_PipeCR)
        call monolis_solver_PipeCR(monoPRM, monoCOM, monoMAT, monoPREC)

      !case (monolis_iter_CABiCGSTAB_noprec)
      !  call monolis_solver_CABiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_PipeBiCGSTAB)
        call monolis_solver_PipeBiCGSTAB(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_PipeBiCGSTAB_noprec)
        call monolis_solver_PipeBiCGSTAB_noprec(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_BiCGSTAB_N128)
        call monolis_solver_BiCGSTAB_N128(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_SOR)
        call monolis_solver_SOR(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_BiCGSAFE)
        call monolis_solver_BiCGSAFE(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_IDRS)
        call monolis_solver_IDRS(monoPRM, monoCOM, monoMAT, monoPREC)

      !case (monolis_iter_GMRES)
      !  call monolis_solver_GMRES(monoPRM, monoCOM, monoMAT)

      case default
        call monolis_std_error_string("monolis_solver_select_R")
        call monolis_std_error_string("please select a linear solver for real numbers")
        call monolis_std_error_stop()
    end select
  end subroutine monolis_solver_select_R

  !> @ingroup dev_solver
  !> 線形ソルバの選択関数
  subroutine monolis_solver_select_C(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_solver_select_R")

    if(monoPRM%Iarray(monolis_prm_I_show_summary) == monolis_I_true .and. monoCOM%my_rank == 0)then
      write(*,"(a)") &
      & "** MONOLIS solver: "//trim(monolis_str_iter(monoPRM%Iarray(monolis_prm_I_method)))//&
      & ", prec: "//trim(monolis_str_prec(monoPRM%Iarray(monolis_prm_I_precond)))
    endif

    select case(monoPRM%Iarray(monolis_prm_I_method))
      case (monolis_iter_COCG)
        call monolis_solver_COCG(monoPRM, monoCOM, monoMAT, monoPREC)

      case default
        call monolis_std_error_string("monolis_solver_select_C")
        call monolis_std_error_string("please select a linear solver for complex numbers")
        call monolis_std_error_stop()
    end select
  end subroutine monolis_solver_select_C
end module mod_monolis_solve
