!> 階層的 VAE (Hierarchical VAE) モジュールのテスト
module mod_monolis_opt_hvae_test
  use mod_monolis
  implicit none

contains

  !> @ingroup optimize
  !> 階層的 VAE 関連サブルーチンの一括テスト
  subroutine monolis_optimize_hvae_test()
    implicit none
    type(monolis_opt_hvae_t) :: hnet
    type(monolis_opt_vae_train_opts) :: opts_local, opts_top
    integer(kint), parameter :: D = 4, H_local = 8, Z_local = 3, H_top = 6, Z_top = 2
    integer(kint), parameter :: N = 32
    real(kdouble_ml) :: X(D, N), zloc(Z_local, N), zloc2(Z_local, N)
    real(kdouble_ml), allocatable :: zall(:,:), ztop(:,:)
    integer(kint) :: i, j, np, N_tot

    !> cover_check 用に対象サブルーチン名を全て登録
    call monolis_std_global_log_string("monolis_opt_hvae_init")
    call monolis_std_global_log_string("monolis_opt_hvae_init_layers")
    call monolis_std_global_log_string("monolis_opt_hvae_finalize")
    call monolis_std_global_log_string("monolis_opt_hvae_gather_latent")
    call monolis_std_global_log_string("monolis_opt_hvae_encode_local")
    call monolis_std_global_log_string("monolis_opt_hvae_encode")
    call monolis_std_global_log_string("monolis_opt_hvae_fit")

    !> 再現性のため Fortran 標準乱数を固定
    call monolis_opt_hvae_test_seed_rng(321)

    np = monolis_mpi_get_global_comm_size()

    !> 初期化と次元チェック (全プロセス共通)
    call monolis_opt_hvae_init(hnet, D, H_local, Z_local, H_top, Z_top)
    call monolis_test_check_eq_I1("hvae_test dim D",       hnet%D,       D)
    call monolis_test_check_eq_I1("hvae_test dim Z_local", hnet%Z_local, Z_local)
    call monolis_test_check_eq_I1("hvae_test dim Z_top",   hnet%Z_top,   Z_top)
    call monolis_test_check_eq_I1("hvae_test local D",     hnet%local%D, D)
    call monolis_test_check_eq_I1("hvae_test local Z",     hnet%local%Z, Z_local)
    if(hnet%my_rank == 0)then
      call monolis_test_check_eq_I1("hvae_test top D", hnet%top%D, Z_local)
      call monolis_test_check_eq_I1("hvae_test top Z", hnet%top%Z, Z_top)
    endif

    !> ダミーデータ (周期的な波形)
    do j = 1, N
      do i = 1, D
        X(i, j) = 0.5d0 + 0.4d0*sin(real(i + j, kdouble))
      end do
    end do

    !> gather_latent 単体テスト: ローカル潜在を既知の値で埋めて集約
    do j = 1, N
      do i = 1, Z_local
        zloc(i, j) = real(hnet%my_rank, kdouble_ml) + 0.01d0*real(i, kdouble_ml)
      end do
    end do
    call monolis_opt_hvae_gather_latent(hnet, zloc, zall)
    call monolis_test_check_eq_I1("hvae_test gather Z dim", size(zall, 1), Z_local)
    if(hnet%my_rank == 0)then
      call monolis_test_check_eq_I1("hvae_test gather N_total", size(zall, 2), np*N)
      !> rank 0 の寄与は先頭ブロックに入る
      call monolis_test_check_eq_R1("hvae_test gather rank0 val", &
        real(zall(1, 1), kdouble), 0.01d0)
    else
      call monolis_test_check_eq_I1("hvae_test gather non-root cols", size(zall, 2), 0)
    endif

    !> 学習 (下位を各プロセスで、上位を rank 0 で)
    opts_local%batch_size = 8
    opts_local%epochs = 5
    opts_local%lr = 1.0e-2_kdouble_ml
    opts_local%r_loss_factor = 1.0_kdouble_ml
    opts_local%kl_warmup_epochs = 2
    opts_local%early_stop_patience = 100
    opts_local%log_every_batches = 0
    opts_local%verbose = .false.

    opts_top%batch_size = 8
    opts_top%epochs = 5
    opts_top%lr = 1.0e-2_kdouble_ml
    opts_top%r_loss_factor = 1.0_kdouble_ml
    opts_top%kl_warmup_epochs = 2
    opts_top%early_stop_patience = 100
    opts_top%log_every_batches = 0
    opts_top%verbose = .false.

    call monolis_opt_hvae_fit(hnet, X, opts_local, opts_top)

    !> 学習済み下位 VAE によるローカルエンコード
    call monolis_opt_hvae_encode_local(hnet, X, zloc2)
    call monolis_test_check_eq_L1("hvae_test encode_local finite", &
      all(zloc2 == zloc2) .and. all(abs(zloc2) < huge(0.0d0)), .true.)

    !> 階層的エンコード (rank 0 で上位潜在を取得)
    call monolis_opt_hvae_encode(hnet, X, ztop)
    call monolis_test_check_eq_I1("hvae_test encode Z_top dim", size(ztop, 1), Z_top)
    if(hnet%my_rank == 0)then
      N_tot = np*N
      call monolis_test_check_eq_I1("hvae_test encode N_total", size(ztop, 2), N_tot)
      call monolis_test_check_eq_L1("hvae_test encode finite", &
        all(ztop == ztop) .and. all(abs(ztop) < huge(0.0d0)), .true.)
    else
      call monolis_test_check_eq_I1("hvae_test encode non-root cols", size(ztop, 2), 0)
    endif

    !> 後始末
    call monolis_opt_hvae_finalize(hnet)
    call monolis_test_check_eq_I1("hvae_test finalize D", hnet%D, 0)

    !> 任意層数の初期化と学習 (下位 enc 2 層 / dec 2 層)
    call monolis_opt_hvae_test_seed_rng(654)
    call monolis_opt_hvae_init_layers(hnet, D, (/ 6, 4 /), Z_local, (/ 5, 7 /), &
      (/ 5 /), Z_top, (/ 5 /))
    call monolis_test_check_eq_I1("hvae_test deep enc layers", size(hnet%local%enc), 2)
    call monolis_test_check_eq_I1("hvae_test deep dec layers", size(hnet%local%dec), 3)
    call monolis_opt_hvae_fit(hnet, X, opts_local, opts_top)
    call monolis_opt_hvae_encode_local(hnet, X, zloc2)
    call monolis_test_check_eq_L1("hvae_test deep encode finite", &
      all(zloc2 == zloc2) .and. all(abs(zloc2) < huge(0.0d0)), .true.)
    call monolis_opt_hvae_finalize(hnet)
  end subroutine monolis_optimize_hvae_test

  !> 再現性確保のための乱数シード
  subroutine monolis_opt_hvae_test_seed_rng(seed)
    implicit none
    integer(kint), intent(in) :: seed
    integer(kint) :: n, i
    integer(kint), allocatable :: s(:)
    call random_seed(size = n)
    allocate(s(n))
    do i = 1, n
      s(i) = seed + 37*i
    end do
    call random_seed(put = s)
    deallocate(s)
  end subroutine monolis_opt_hvae_test_seed_rng

end module mod_monolis_opt_hvae_test
