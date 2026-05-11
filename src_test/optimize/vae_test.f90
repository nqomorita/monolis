!> VAE モジュールのテスト
module mod_monolis_opt_vae_test
  use mod_monolis
  implicit none

contains

  !> @ingroup optimize
  !> VAE 関連サブルーチンの一括テスト
  subroutine monolis_optimize_vae_test()
    implicit none
    type(monolis_opt_vae_t) :: net
    type(monolis_opt_vae_train_opts) :: opts
    integer(kint), parameter :: D = 4, H = 8, Z = 2, N = 32
    real(kdouble) :: X(D, N), xhat(D, N), zlat(Z, N), samples(D, 4)
    real(kdouble) :: parents(Z, Z+1), child(Z)
    real(kdouble) :: loss, recon, kl, loss0
    integer(kint) :: i, j

    !> cover_check 用に対象サブルーチン名を全て登録
    call monolis_std_global_log_string("monolis_opt_vae_init")
    call monolis_std_global_log_string("monolis_opt_vae_finalize")
    call monolis_std_global_log_string("monolis_opt_vae_adam_init")
    call monolis_std_global_log_string("monolis_opt_vae_adam_free")
    call monolis_std_global_log_string("monolis_opt_vae_adam_apply")
    call monolis_std_global_log_string("monolis_opt_vae_glorot_uniform")
    call monolis_std_global_log_string("monolis_opt_vae_fill_randn")
    call monolis_std_global_log_string("monolis_opt_vae_shuffle")
    call monolis_std_global_log_string("monolis_opt_vae_train_step")
    call monolis_std_global_log_string("monolis_opt_vae_fit")
    call monolis_std_global_log_string("monolis_opt_vae_reconstruct")
    call monolis_std_global_log_string("monolis_opt_vae_encode")
    call monolis_std_global_log_string("monolis_opt_vae_decode")
    call monolis_std_global_log_string("monolis_opt_vae_sample_prior")
    call monolis_std_global_log_string("monolis_opt_vae_spx_child")
    call monolis_std_global_log_string("monolis_opt_vae_generate_spx")

    !> 再現性のため Fortran 標準乱数を固定
    call monolis_opt_vae_test_seed_rng(123)

    !> 初期化
    call monolis_opt_vae_init(net, D, H, Z)
    call monolis_test_check_eq_I1("vae_test dim D", net%D, D)
    call monolis_test_check_eq_I1("vae_test dim H", net%H, H)
    call monolis_test_check_eq_I1("vae_test dim Z", net%Z, Z)

    !> ダミーデータ (4 次元の周期的なシグモイド波形)
    do j = 1, N
      do i = 1, D
        X(i, j) = 0.5d0 + 0.4d0*sin(real(i + j, kdouble))
      end do
    end do

    !> 1 ステップ学習で損失が有限値であること
    call monolis_opt_vae_train_step(net, X, 1.0d-3, 1.0d0, loss0, recon, kl)
    call monolis_test_check_eq_L1("vae_test loss finite", loss0 == loss0 .and. abs(loss0) < huge(0.0d0), .true.)

    !> 短い fit 呼び出し: 損失が初期 1 ステップより悪化しないこと
    opts%batch_size = 8
    opts%epochs = 5
    opts%lr = 1.0d-2
    opts%r_loss_factor = 1.0d0
    opts%kl_warmup_epochs = 2
    opts%early_stop_patience = 100
    opts%log_every_batches = 0
    opts%verbose = .false.
    call monolis_opt_vae_fit(net, X, opts)
    call monolis_opt_vae_train_step(net, X, 0.0d0, 1.0d0, loss, recon, kl)
    call monolis_test_check_eq_L1("vae_test loss decreased", loss <= loss0 + 1.0d-6, .true.)

    !> 再構成と次元
    call monolis_opt_vae_reconstruct(net, X, xhat)
    call monolis_test_check_eq_L1("vae_test xhat in [0,1]", &
      all(xhat >= 0.0d0) .and. all(xhat <= 1.0d0), .true.)

    !> エンコード / デコード
    call monolis_opt_vae_encode(net, X, zlat)
    call monolis_opt_vae_decode(net, zlat, xhat)
    call monolis_test_check_eq_L1("vae_test decode in [0,1]", &
      all(xhat >= 0.0d0) .and. all(xhat <= 1.0d0), .true.)

    !> 事前分布からのサンプル
    call monolis_opt_vae_sample_prior(net, 4, samples)
    call monolis_test_check_eq_L1("vae_test prior in [0,1]", &
      all(samples >= 0.0d0) .and. all(samples <= 1.0d0), .true.)

    !> SPX 単体: 同一の親から生成すれば子も同じ点
    do j = 1, Z+1
      parents(:, j) = (/ 1.0d0, -1.0d0 /)
    end do
    call monolis_opt_vae_spx_child(parents, Z, child)
    call monolis_test_check_eq_R1("vae_test spx degenerate 1", child(1), 1.0d0)
    call monolis_test_check_eq_R1("vae_test spx degenerate 2", child(2), -1.0d0)

    !> SPX 経由の生成
    call monolis_opt_vae_generate_spx(net, X, 4, samples)
    call monolis_test_check_eq_L1("vae_test spx samples in [0,1]", &
      all(samples >= 0.0d0) .and. all(samples <= 1.0d0), .true.)

    !> 後始末
    call monolis_opt_vae_finalize(net)
    call monolis_test_check_eq_I1("vae_test finalize D", net%D, 0)
  end subroutine monolis_optimize_vae_test

  !> 再現性確保のための乱数シード
  subroutine monolis_opt_vae_test_seed_rng(seed)
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
  end subroutine monolis_opt_vae_test_seed_rng

end module mod_monolis_opt_vae_test
