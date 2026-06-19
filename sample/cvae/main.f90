!> CVAE (Conditional VAE) サンプル
!>
!> 概要:
!>   2 つの 1 次元周期パターン (sin/cos) からなる合成データを、
!>   クラスラベルを条件ベクトル (one-hot) として与えて学習する。
!>   学習後は条件付き再構成と、条件を指定した事前分布サンプリング
!>   (条件付き生成) を実行し、再構成 RMSE を標準出力に表示する。
!>
!> 使い方:
!>   ./run.sh
program main
  use mod_monolis
  use mod_monolis_utils
  use mod_monolis_def_opt
  implicit none

  integer(kint), parameter :: D = 16, C = 2, H = 32, Z = 4, N = 256
  real(kdouble), parameter :: PI = 3.141592653589793d0
  type(monolis_opt_cvae_t) :: net
  type(monolis_opt_vae_train_opts) :: opts
  real(kdouble_ml) :: X(D, N), Cnd(C, N), xhat(D, N)
  real(kdouble_ml) :: Cgen(C, 8), gen(D, 8)
  real(kdouble_ml) :: rmse
  integer(kint) :: i, j, label, seed_size
  integer(kint), allocatable :: seed_arr(:)
  real(kdouble) :: r, val

  call monolis_global_initialize()

  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,'(a)') "================================================"
    write(*,'(a)') " CVAE sample : conditioned 1D wave dataset"
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') "  D=", D, " C=", C, " H=", H, " Z=", Z, " N=", N
    write(*,'(a)') "================================================"
  endif

  !> 乱数シード固定
  call random_seed(size = seed_size)
  allocate(seed_arr(seed_size))
  do i = 1, seed_size
    seed_arr(i) = 2025 + 13*i
  end do
  call random_seed(put = seed_arr)
  deallocate(seed_arr)

  !> データ生成: ラベル 0 (sin 波) とラベル 1 (cos 波) をランダムに混在
  !> 条件 Cnd はラベルの one-hot ベクトル (C=2)
  do j = 1, N
    call random_number(r)
    label = 0
    if(r > 0.5d0) label = 1
    Cnd(:, j) = 0.0_kdouble_ml
    Cnd(label + 1, j) = 1.0_kdouble_ml
    do i = 1, D
      if(label == 0)then
        val = 0.5d0 + 0.4d0 * sin(2.0d0*PI*real(i, kdouble)/real(D, kdouble))
      else
        val = 0.5d0 + 0.4d0 * cos(2.0d0*PI*real(i, kdouble)/real(D, kdouble))
      endif
      X(i, j) = real(val, kdouble_ml)
    end do
  end do

  !> モデル構築と学習
  call monolis_opt_cvae_init(net, D, C, H, Z)

  opts%batch_size = 32
  opts%epochs = 50
  opts%lr = 5.0d-3
  opts%r_loss_factor = 1.0d0
  opts%kl_warmup_epochs = 10
  opts%early_stop_patience = 30
  opts%log_every_batches = 0
  opts%verbose = .true.

  call monolis_opt_cvae_fit(net, X, Cnd, opts)

  !> 条件付き再構成
  call monolis_opt_cvae_reconstruct(net, X, Cnd, xhat)
  rmse = sqrt(sum((X - xhat)**2) / real(D*N, kdouble_ml))
  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,'(a,es12.4)') " reconstruction RMSE = ", rmse
  endif

  !> 条件付き生成: ラベル 0 (sin) を 4 個、ラベル 1 (cos) を 4 個
  do j = 1, 8
    Cgen(:, j) = 0.0_kdouble_ml
    if(j <= 4)then
      Cgen(1, j) = 1.0_kdouble_ml
    else
      Cgen(2, j) = 1.0_kdouble_ml
    endif
  end do
  call monolis_opt_cvae_sample_prior(net, Cgen, gen)
  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,'(a)') " --- conditioned prior samples (first 4 dims of each sample) ---"
    do j = 1, 8
      if(j <= 4)then
        write(*,'(a,i0,a,4f8.4)') "  label 0 (sin) sample ", j, " : ", gen(1:4, j)
      else
        write(*,'(a,i0,a,4f8.4)') "  label 1 (cos) sample ", j-4, " : ", gen(1:4, j)
      endif
    end do
  endif

  call monolis_opt_cvae_finalize(net)
  call monolis_global_finalize()
end program main
