!> 階層的変分オートエンコーダ (Hierarchical VAE) ライブラリ
!> @details MPI 並列で動作する 2 階層の VAE。
!>          1. 各 MPI プロセスがローカルデータに対して通常の VAE
!>             (mod_monolis_opt_vae) を学習する (下位層)。
!>          2. 学習済みローカル VAE で入力を潜在ベクトル (latent) に変換する。
!>          3. 全プロセスの latent を MPI rank 0 へ全集約する。
!>          4. rank 0 のみが集約 latent を入力データとして上位層 VAE を学習する。
!>          下位層・上位層ともに mod_monolis_opt_vae の実装を再利用する。
!>          機械学習に関わる実数は 32bit 浮動小数点 (kdouble_ml) で計算する。
module mod_monolis_opt_hvae
  use mod_monolis_utils
  use mod_monolis_def_opt
  use mod_monolis_opt_vae_util
  use mod_monolis_opt_vae
  implicit none

  private
  public :: monolis_opt_hvae_t
  public :: monolis_opt_hvae_init
  public :: monolis_opt_hvae_init_layers
  public :: monolis_opt_hvae_finalize
  public :: monolis_opt_hvae_gather_latent
  public :: monolis_opt_hvae_encode_local
  public :: monolis_opt_hvae_encode
  public :: monolis_opt_hvae_fit

  !> @ingroup optimize
  !> 階層的 VAE モデル (下位ローカル VAE + 上位 VAE)
  type :: monolis_opt_hvae_t
    !> [in] 入力次元数 (下位 VAE の入力次元)
    integer(kint) :: D = 0
    !> [in] 下位 VAE の潜在次元数 (= 上位 VAE の入力次元)
    integer(kint) :: Z_local = 0
    !> [in] 上位 VAE の潜在次元数
    integer(kint) :: Z_top = 0
    !> [in] MPI コミュニケータ
    integer(kint) :: comm = 0
    !> [in] 自プロセスの MPI ランク
    integer(kint) :: my_rank = 0
    !> [in] MPI プロセス数
    integer(kint) :: comm_size = 1
    !> [in,out] 下位 VAE (全プロセスが保持)
    type(monolis_opt_vae_t) :: local
    !> [in,out] 上位 VAE (rank 0 のみ初期化・学習する)
    type(monolis_opt_vae_t) :: top
  end type monolis_opt_hvae_t

contains

  !> @ingroup optimize
  !> 階層的 VAE を初期化する (旧 API。各層 1 隠れ層で構成)
  subroutine monolis_opt_hvae_init(hnet, D, H_local, Z_local, H_top, Z_top)
    implicit none
    !> [out] 初期化対象の階層的 VAE モデル
    type(monolis_opt_hvae_t), intent(out) :: hnet
    !> [in] 入力次元数
    integer(kint), intent(in) :: D
    !> [in] 下位 VAE の隠れ層次元数
    integer(kint), intent(in) :: H_local
    !> [in] 下位 VAE の潜在次元数
    integer(kint), intent(in) :: Z_local
    !> [in] 上位 VAE の隠れ層次元数
    integer(kint), intent(in) :: H_top
    !> [in] 上位 VAE の潜在次元数
    integer(kint), intent(in) :: Z_top

    call monolis_opt_hvae_init_layers(hnet, D, (/ H_local /), Z_local, (/ H_local /), &
      (/ H_top /), Z_top, (/ H_top /))
  end subroutine monolis_opt_hvae_init

  !> @ingroup optimize
  !> 任意層数の階層的 VAE を初期化する
  !> @details 下位 VAE は全プロセスで初期化する。上位 VAE は rank 0 のみ初期化する。
  !>          上位 VAE の入力次元は下位 VAE の潜在次元 Z_local に一致する。
  subroutine monolis_opt_hvae_init_layers(hnet, D, H_enc_local, Z_local, H_dec_local, &
      H_enc_top, Z_top, H_dec_top)
    implicit none
    !> [out] 初期化対象の階層的 VAE モデル
    type(monolis_opt_hvae_t), intent(out) :: hnet
    !> [in] 入力次元数
    integer(kint), intent(in) :: D
    !> [in] 下位 VAE のエンコーダ隠れ層次元の列 (size>=1)
    integer(kint), intent(in) :: H_enc_local(:)
    !> [in] 下位 VAE の潜在次元数
    integer(kint), intent(in) :: Z_local
    !> [in] 下位 VAE のデコーダ隠れ層次元の列 (size>=1)
    integer(kint), intent(in) :: H_dec_local(:)
    !> [in] 上位 VAE のエンコーダ隠れ層次元の列 (size>=1)
    integer(kint), intent(in) :: H_enc_top(:)
    !> [in] 上位 VAE の潜在次元数
    integer(kint), intent(in) :: Z_top
    !> [in] 上位 VAE のデコーダ隠れ層次元の列 (size>=1)
    integer(kint), intent(in) :: H_dec_top(:)

    call monolis_std_debug_log_header("monolis_opt_hvae_init_layers")

    hnet%D         = D
    hnet%Z_local   = Z_local
    hnet%Z_top     = Z_top
    hnet%comm      = monolis_mpi_get_global_comm()
    hnet%my_rank   = monolis_mpi_get_global_my_rank()
    hnet%comm_size = monolis_mpi_get_global_comm_size()

    !> 下位 VAE は全プロセスで初期化
    call monolis_opt_vae_init_layers(hnet%local, D, H_enc_local, Z_local, H_dec_local)

    !> 上位 VAE は rank 0 のみ初期化 (入力次元は下位潜在次元 Z_local)
    if(hnet%my_rank == 0)then
      call monolis_opt_vae_init_layers(hnet%top, Z_local, H_enc_top, Z_top, H_dec_top)
    endif
  end subroutine monolis_opt_hvae_init_layers

  !> @ingroup optimize
  !> 階層的 VAE が保持するメモリを解放する
  subroutine monolis_opt_hvae_finalize(hnet)
    implicit none
    !> [in,out] 解放対象の階層的 VAE モデル
    type(monolis_opt_hvae_t), intent(inout) :: hnet

    call monolis_opt_vae_finalize(hnet%local)
    if(hnet%my_rank == 0)then
      call monolis_opt_vae_finalize(hnet%top)
    endif
    hnet%D = 0; hnet%Z_local = 0; hnet%Z_top = 0
    hnet%comm = 0; hnet%my_rank = 0; hnet%comm_size = 1
  end subroutine monolis_opt_hvae_finalize

  !> @ingroup optimize
  !> 各プロセスのローカル潜在ベクトルを rank 0 へ全集約する
  !> @details 各プロセスの zloc (Z x N_local) を rank 0 に集め、
  !>          rank 0 では zall (Z x N_total) を返す (列はランク順に連結)。
  !>          rank 0 以外では zall は (Z x 0) で返る。
  subroutine monolis_opt_hvae_gather_latent(hnet, zloc, zall)
    implicit none
    !> [in] 階層的 VAE モデル
    type(monolis_opt_hvae_t), intent(in) :: hnet
    !> [in] ローカル潜在ベクトル (Z x N_local)
    real(kdouble_ml), intent(in) :: zloc(:,:)
    !> [out] 全集約された潜在ベクトル (rank 0: Z x N_total、他: Z x 0)
    real(kdouble_ml), allocatable, intent(out) :: zall(:,:)
    integer(kint) :: Z, N_loc, sc, np, total, N_tot, i
    integer(kint), allocatable :: rc(:), disp(:)
    real(kdouble_ml), allocatable :: sbuf(:), rbuf(:)

    call monolis_std_debug_log_header("monolis_opt_hvae_gather_latent")

    Z     = size(zloc, 1)
    N_loc = size(zloc, 2)
    np    = hnet%comm_size
    sc    = Z * N_loc

    !> 各プロセスの送信個数を全プロセスで共有し、格納位置 (disp) を計算
    call monolis_alloc_I_1d(rc, np)
    call monolis_alloc_I_1d(disp, np)
    call monolis_allgather_I1(sc, rc, hnet%comm)
    disp(1) = 0
    do i = 2, np
      disp(i) = disp(i-1) + rc(i-1)
    end do
    total = disp(np) + rc(np)

    !> 送信バッファ (列優先で平坦化)
    call monolis_alloc_F_1d(sbuf, max(sc, 1))
    if(sc > 0) sbuf(1:sc) = reshape(zloc, (/ sc /))

    !> 受信バッファ (rank 0 のみ全長を確保)
    if(hnet%my_rank == 0)then
      call monolis_alloc_F_1d(rbuf, max(total, 1))
    else
      call monolis_alloc_F_1d(rbuf, 1)
    endif

    call monolis_gather_V_F(sbuf, sc, rbuf, rc, disp, 0, hnet%comm)

    if(hnet%my_rank == 0)then
      N_tot = total / Z
      call monolis_alloc_F_2d(zall, Z, N_tot)
      if(total > 0) zall = reshape(rbuf(1:total), (/ Z, N_tot /))
    else
      call monolis_alloc_F_2d(zall, Z, 0)
    endif
  end subroutine monolis_opt_hvae_gather_latent

  !> @ingroup optimize
  !> 学習済みローカル VAE で入力をローカル潜在ベクトルにエンコードする
  subroutine monolis_opt_hvae_encode_local(hnet, X, zloc)
    implicit none
    !> [in] 階層的 VAE モデル
    type(monolis_opt_hvae_t), intent(in) :: hnet
    !> [in] ローカル入力データ (D x N_local)
    real(kdouble_ml), intent(in) :: X(:,:)
    !> [out] ローカル潜在ベクトル (Z_local x N_local)
    real(kdouble_ml), intent(out) :: zloc(:,:)

    call monolis_std_debug_log_header("monolis_opt_hvae_encode_local")
    call monolis_opt_vae_encode(hnet%local, X, zloc)
  end subroutine monolis_opt_hvae_encode_local

  !> @ingroup optimize
  !> 入力を階層的にエンコードして上位潜在ベクトルを得る
  !> @details ローカル VAE でローカル潜在に変換 → rank 0 へ全集約 →
  !>          rank 0 で上位 VAE エンコード。上位潜在 ztop は rank 0 のみ
  !>          (Z_top x N_total)、他プロセスでは (Z_top x 0) で返る。
  subroutine monolis_opt_hvae_encode(hnet, X, ztop)
    implicit none
    !> [in] 階層的 VAE モデル
    type(monolis_opt_hvae_t), intent(in) :: hnet
    !> [in] ローカル入力データ (D x N_local)
    real(kdouble_ml), intent(in) :: X(:,:)
    !> [out] 上位潜在ベクトル (rank 0: Z_top x N_total、他: Z_top x 0)
    real(kdouble_ml), allocatable, intent(out) :: ztop(:,:)
    real(kdouble_ml), allocatable :: zloc(:,:), zall(:,:)
    integer(kint) :: N_tot

    call monolis_std_debug_log_header("monolis_opt_hvae_encode")

    !> ローカルエンコード → 全集約
    call monolis_alloc_F_2d(zloc, hnet%Z_local, size(X, 2))
    call monolis_opt_hvae_encode_local(hnet, X, zloc)
    call monolis_opt_hvae_gather_latent(hnet, zloc, zall)

    !> rank 0 で上位エンコード
    if(hnet%my_rank == 0)then
      N_tot = size(zall, 2)
      call monolis_alloc_F_2d(ztop, hnet%Z_top, N_tot)
      call monolis_opt_vae_encode(hnet%top, zall, ztop)
    else
      call monolis_alloc_F_2d(ztop, hnet%Z_top, 0)
    endif
  end subroutine monolis_opt_hvae_encode

  !> @ingroup optimize
  !> 階層的 VAE の学習
  !> @details 手順:
  !>          1. 各プロセスがローカルデータ X で下位 VAE を学習する。
  !>          2. 学習済み下位 VAE で X をローカル潜在に変換する。
  !>          3. ローカル潜在を rank 0 へ全集約する。
  !>          4. rank 0 のみが集約潜在を入力として上位 VAE を学習する。
  subroutine monolis_opt_hvae_fit(hnet, X, opts_local, opts_top)
    implicit none
    !> [in,out] 学習対象の階層的 VAE モデル
    type(monolis_opt_hvae_t), intent(inout) :: hnet
    !> [in] ローカル学習データ (D x N_local)
    real(kdouble_ml), intent(in) :: X(:,:)
    !> [in] 下位 VAE の学習オプション
    type(monolis_opt_vae_train_opts), intent(in) :: opts_local
    !> [in] 上位 VAE の学習オプション
    type(monolis_opt_vae_train_opts), intent(in) :: opts_top
    real(kdouble_ml), allocatable :: zloc(:,:), zall(:,:)

    call monolis_std_debug_log_header("monolis_opt_hvae_fit")

    !> 1. 下位 (ローカル) VAE を各プロセスで学習
    call monolis_opt_vae_fit(hnet%local, X, opts_local)

    !> 2. 入力をローカル潜在に変換
    call monolis_alloc_F_2d(zloc, hnet%Z_local, size(X, 2))
    call monolis_opt_hvae_encode_local(hnet, X, zloc)

    !> 3. ローカル潜在を rank 0 へ全集約
    call monolis_opt_hvae_gather_latent(hnet, zloc, zall)

    !> 4. rank 0 のみで上位 VAE を学習
    if(hnet%my_rank == 0)then
      call monolis_opt_vae_fit(hnet%top, zall, opts_top)
    endif
  end subroutine monolis_opt_hvae_fit

end module mod_monolis_opt_hvae
