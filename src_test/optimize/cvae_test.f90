!> Conditional VAE モジュールのテスト
module mod_monolis_opt_cvae_test
  use mod_monolis
  implicit none

contains

  !> @ingroup optimize
  !> Conditional VAE 関連サブルーチンの一括テスト
  subroutine monolis_optimize_cvae_test()
    implicit none
    type(monolis_opt_cvae_t) :: net
    type(monolis_opt_vae_train_opts) :: opts
    integer(kint), parameter :: D = 4, C = 2, H = 8, Z = 2, N = 32
    real(kdouble_ml) :: X(D, N), Cnd(C, N), xhat(D, N), zlat(Z, N)
    real(kdouble_ml) :: Cgen(C, 4), samples(D, 4)
    real(kdouble_ml) :: loss, recon, kl, loss0
    integer(kint) :: i, j

    !> cover_check 用に対象サブルーチン名を全て登録
    call monolis_std_global_log_string("monolis_opt_cvae_concat_rows")
    call monolis_std_global_log_string("monolis_opt_cvae_init")
    call monolis_std_global_log_string("monolis_opt_cvae_init_layers")
    call monolis_std_global_log_string("monolis_opt_cvae_finalize")
    call monolis_std_global_log_string("monolis_opt_cvae_train_step")
    call monolis_std_global_log_string("monolis_opt_cvae_fit")
    call monolis_std_global_log_string("monolis_opt_cvae_encode_mu_lv")
    call monolis_std_global_log_string("monolis_opt_cvae_decode_alloc")
    call monolis_std_global_log_string("monolis_opt_cvae_reconstruct")
    call monolis_std_global_log_string("monolis_opt_cvae_encode")
    call monolis_std_global_log_string("monolis_opt_cvae_decode")
    call monolis_std_global_log_string("monolis_opt_cvae_sample_prior")

    !> 再現性のため Fortran 標準乱数を固定
    call monolis_opt_cvae_test_seed_rng(321)

    !> 初期化
    call monolis_opt_cvae_init(net, D, C, H, Z)
    call monolis_test_check_eq_I1("cvae_test dim D", net%D, D)
    call monolis_test_check_eq_I1("cvae_test dim C", net%C, C)
    call monolis_test_check_eq_I1("cvae_test dim H", net%H, H)
    call monolis_test_check_eq_I1("cvae_test dim Z", net%Z, Z)
    call monolis_test_check_eq_I1("cvae_test enc(1) in", net%enc(1)%in_dim, D + C)
    call monolis_test_check_eq_I1("cvae_test dec(1) in", net%dec(1)%in_dim, Z + C)

    !> ダミーデータ (D 次元の周期波形) と 2 次元の条件
    do j = 1, N
      do i = 1, D
        X(i, j) = 0.5d0 + 0.4d0*sin(real(i + j, kdouble))
      end do
      Cnd(1, j) = 0.5d0 + 0.5d0*sin(real(j, kdouble))
      Cnd(2, j) = 0.5d0 + 0.5d0*cos(real(j, kdouble))
    end do

    !> 1 ステップ学習で損失が有限値であること
    call monolis_opt_cvae_train_step(net, X, Cnd, 1.0e-3_kdouble_ml, 1.0_kdouble_ml, loss0, recon, kl)
    call monolis_test_check_eq_L1("cvae_test loss finite", &
      loss0 == loss0 .and. abs(loss0) < huge(0.0d0), .true.)

    !> 短い fit 呼び出し: 損失が初期 1 ステップより悪化しないこと
    opts%batch_size = 8
    opts%epochs = 5
    opts%lr = 1.0e-2_kdouble_ml
    opts%r_loss_factor = 1.0_kdouble_ml
    opts%kl_warmup_epochs = 2
    opts%early_stop_patience = 100
    opts%log_every_batches = 0
    opts%verbose = .false.
    call monolis_opt_cvae_fit(net, X, Cnd, opts)
    call monolis_opt_cvae_train_step(net, X, Cnd, 0.0_kdouble_ml, 1.0_kdouble_ml, loss, recon, kl)
    call monolis_test_check_eq_L1("cvae_test loss decreased", loss <= loss0 + 1.0d-6, .true.)

    !> 再構成と次元
    call monolis_opt_cvae_reconstruct(net, X, Cnd, xhat)
    call monolis_test_check_eq_L1("cvae_test xhat in [0,1]", &
      all(xhat >= 0.0d0) .and. all(xhat <= 1.0d0), .true.)

    !> エンコード / デコード
    call monolis_opt_cvae_encode(net, X, Cnd, zlat)
    call monolis_opt_cvae_decode(net, zlat, Cnd, xhat)
    call monolis_test_check_eq_L1("cvae_test decode in [0,1]", &
      all(xhat >= 0.0d0) .and. all(xhat <= 1.0d0), .true.)

    !> 条件付き事前分布からのサンプル
    do j = 1, 4
      Cgen(1, j) = 0.25d0*real(j, kdouble)
      Cgen(2, j) = 1.0d0 - 0.25d0*real(j, kdouble)
    end do
    call monolis_opt_cvae_sample_prior(net, Cgen, samples)
    call monolis_test_check_eq_L1("cvae_test prior in [0,1]", &
      all(samples >= 0.0d0) .and. all(samples <= 1.0d0), .true.)

    !> 後始末
    call monolis_opt_cvae_finalize(net)
    call monolis_test_check_eq_I1("cvae_test finalize D", net%D, 0)
    call monolis_test_check_eq_I1("cvae_test finalize C", net%C, 0)

    !> 任意層数 (隠れ 2 層 enc / 2 層 dec) の初期化と 1 ステップ学習
    call monolis_opt_cvae_test_seed_rng(654)
    call monolis_opt_cvae_init_layers(net, D, C, (/ 6, 4 /), Z, (/ 5, 7 /))
    call monolis_test_check_eq_I1("cvae_test deep enc layers", size(net%enc), 2)
    call monolis_test_check_eq_I1("cvae_test deep dec layers", size(net%dec), 3)
    call monolis_test_check_eq_I1("cvae_test deep enc(1) in", net%enc(1)%in_dim, D + C)
    call monolis_test_check_eq_I1("cvae_test deep dec(1) in", net%dec(1)%in_dim, Z + C)
    call monolis_test_check_eq_I1("cvae_test deep dec(last) out", net%dec(3)%out_dim, D)
    call monolis_opt_cvae_train_step(net, X, Cnd, 1.0e-3_kdouble_ml, 1.0_kdouble_ml, loss, recon, kl)
    call monolis_test_check_eq_L1("cvae_test deep loss finite", &
      loss == loss .and. abs(loss) < huge(0.0d0), .true.)
    call monolis_opt_cvae_reconstruct(net, X, Cnd, xhat)
    call monolis_test_check_eq_L1("cvae_test deep xhat in [0,1]", &
      all(xhat >= 0.0d0) .and. all(xhat <= 1.0d0), .true.)
    call monolis_opt_cvae_finalize(net)
  end subroutine monolis_optimize_cvae_test

  !> 再現性確保のための乱数シード
  subroutine monolis_opt_cvae_test_seed_rng(seed)
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
  end subroutine monolis_opt_cvae_test_seed_rng

end module mod_monolis_opt_cvae_test
