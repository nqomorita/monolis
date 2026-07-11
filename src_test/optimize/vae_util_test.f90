!> VAE / CVAE 共通ユーティリティモジュールのテスト
module mod_monolis_opt_vae_util_test
  use mod_monolis
  implicit none

contains

  !> @ingroup optimize
  !> vae_util 関連サブルーチンの一括テスト
  subroutine monolis_optimize_vae_util_test()
    implicit none
    type(monolis_opt_vae_layer_t) :: layers(2)
    type(monolis_opt_vae_cache_t), allocatable :: cache(:)
    type(monolis_opt_vae_grad_t),  allocatable :: grads(:)
    real(kdouble_ml) :: X(3, 5), dY(2, 5), randn(4, 4)
    real(kdouble_ml) :: parents(2, 3), child(2)
    real(kdouble_ml), allocatable :: top(:,:), dX_in(:,:)
    integer(kint) :: perm(6), i, j, cnt

    !> cover_check 用に対象サブルーチン名を全て登録
    call monolis_std_global_log_string("monolis_opt_vae_layer_init")
    call monolis_std_global_log_string("monolis_opt_vae_layer_free")
    call monolis_std_global_log_string("monolis_opt_vae_glorot_uniform")
    call monolis_std_global_log_string("monolis_opt_vae_fill_randn")
    call monolis_std_global_log_string("monolis_opt_vae_shuffle")
    call monolis_std_global_log_string("monolis_opt_vae_layer_forward")
    call monolis_std_global_log_string("monolis_opt_vae_layer_forward_kernel")
    call monolis_std_global_log_string("monolis_opt_vae_forward_stack")
    call monolis_std_global_log_string("monolis_opt_vae_top_post")
    call monolis_std_global_log_string("monolis_opt_vae_layer_backward")
    call monolis_std_global_log_string("monolis_opt_vae_layer_backward_kernel")
    call monolis_std_global_log_string("monolis_opt_vae_backward_stack")
    call monolis_std_global_log_string("monolis_opt_vae_dev_copy2")
    call monolis_std_global_log_string("monolis_opt_vae_cache_free")
    call monolis_std_global_log_string("monolis_opt_vae_grads_free")
    call monolis_std_global_log_string("monolis_opt_vae_spx_child")

    !> 入力データ
    do j = 1, 5
      do i = 1, 3
        X(i, j) = 0.1d0*real(i, kdouble) - 0.2d0*real(j, kdouble)
      end do
    end do

    !> 2 層 (3 -> 4 ReLU -> 2 sigmoid) の順伝播・逆伝播
    call monolis_opt_vae_layer_init(layers(1), 3, 4, monolis_opt_vae_act_relu)
    call monolis_opt_vae_layer_init(layers(2), 4, 2, monolis_opt_vae_act_sigmoid)
    call monolis_test_check_eq_I1("vae_util in_dim", layers(1)%in_dim, 3)
    call monolis_test_check_eq_I1("vae_util out_dim", layers(2)%out_dim, 2)

    call monolis_opt_vae_forward_stack(layers, X, cache)
    call monolis_test_check_eq_I1("vae_util cache size", size(cache), 2)
    top = monolis_opt_vae_top_post(cache)
    call monolis_test_check_eq_I1("vae_util top rows", size(top, 1), 2)
    call monolis_test_check_eq_L1("vae_util top in [0,1]", &
      all(top >= 0.0d0) .and. all(top <= 1.0d0), .true.)

    dY = 1.0d0
    call monolis_opt_vae_backward_stack(layers, X, cache, dY, grads, dX_in)
    call monolis_test_check_eq_I1("vae_util grads size", size(grads), 2)
    call monolis_test_check_eq_I1("vae_util gW1 rows", size(grads(1)%gW, 1), 3)
    call monolis_test_check_eq_I1("vae_util dX rows", size(dX_in, 1), 3)
    call monolis_test_check_eq_L1("vae_util grad finite", &
      all(grads(1)%gW == grads(1)%gW), .true.)

    call monolis_opt_vae_cache_free(cache)
    call monolis_opt_vae_grads_free(grads)
    call monolis_test_check_eq_L1("vae_util cache freed", allocated(cache), .false.)
    call monolis_test_check_eq_L1("vae_util grads freed", allocated(grads), .false.)
    call monolis_opt_vae_layer_free(layers(1))
    call monolis_opt_vae_layer_free(layers(2))
    call monolis_test_check_eq_I1("vae_util layer freed", layers(1)%in_dim, 0)

    !> 標準正規乱数
    call monolis_opt_vae_fill_randn(randn)
    call monolis_test_check_eq_L1("vae_util randn finite", &
      all(randn == randn) .and. all(abs(randn) < huge(0.0d0)), .true.)

    !> シャッフルが順列を保つこと
    do i = 1, 6
      perm(i) = i
    end do
    call monolis_opt_vae_shuffle(perm)
    cnt = 0
    do i = 1, 6
      do j = 1, 6
        if(perm(j) == i) cnt = cnt + 1
      end do
    end do
    call monolis_test_check_eq_I1("vae_util shuffle permutation", cnt, 6)

    !> SPX 単体: 同一の親から生成すれば子も同じ点
    do j = 1, 3
      parents(:, j) = (/ 2.0d0, -3.0d0 /)
    end do
    call monolis_opt_vae_spx_child(parents, 2, child)
    call monolis_test_check_eq_R1("vae_util spx degenerate 1", real(child(1), kdouble), 2.0d0)
    call monolis_test_check_eq_R1("vae_util spx degenerate 2", real(child(2), kdouble), -3.0d0)
  end subroutine monolis_optimize_vae_util_test

end module mod_monolis_opt_vae_util_test
