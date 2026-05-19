!> VAE サンプル
!>
!> 概要:
!>   2 つの 1 次元周期パターン (sin/cos) からなる合成データを学習し、
!>   再構成・潜在エンコード・事前分布サンプリング・SPX 交叉サンプリング
!>   を実行する。結果は標準出力に再構成 RMSE を表示する。
!>
!> 使い方:
!>   ./run.sh
program main
  use mod_monolis
  use mod_monolis_utils
  implicit none

  integer(kint), parameter :: D = 16, H = 32, Z = 4, N = 256
  type(monolis_opt_vae_t) :: net
  type(monolis_opt_vae_train_opts) :: opts
  real(kdouble) :: X(D, N), xhat(D, N), samples(D, 8)
  real(kdouble) :: rmse
  integer(kint) :: i, j, label, seed_size
  integer(kint), allocatable :: seed_arr(:)
  real(kdouble) :: r

  call monolis_global_initialize()

  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,'(a)') "================================================"
    write(*,'(a)') " VAE sample : synthetic 1D wave dataset"
    write(*,'(a,i0,a,i0,a,i0,a,i0)') "  D=", D, " H=", H, " Z=", Z, " N=", N
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
  do j = 1, N
    call random_number(r)
    label = 0
    if(r > 0.5d0) label = 1
    do i = 1, D
      if(label == 0)then
        X(i, j) = 0.5d0 + 0.4d0 * sin(2.0d0*3.141592653589793d0*real(i, kdouble)/real(D, kdouble))
      else
        X(i, j) = 0.5d0 + 0.4d0 * cos(2.0d0*3.141592653589793d0*real(i, kdouble)/real(D, kdouble))
      endif
    end do
  end do

  !> モデル構築と学習
  call monolis_opt_vae_init(net, D, H, Z)

  opts%batch_size = 32
  opts%epochs = 50
  opts%lr = 5.0d-3
  opts%r_loss_factor = 1.0d0
  opts%kl_warmup_epochs = 10
  opts%early_stop_patience = 30
  opts%log_every_batches = 0
  opts%verbose = .true.

  call monolis_opt_vae_fit(net, X, opts)

  !> 再構成
  call monolis_opt_vae_reconstruct(net, X, xhat)
  rmse = sqrt(sum((X - xhat)**2) / real(D*N, kdouble))
  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,'(a,es12.4)') " reconstruction RMSE = ", rmse
  endif

  !> 事前分布から 8 個サンプリング
  call monolis_opt_vae_sample_prior(net, 8, samples)
  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,'(a)') " --- prior samples (first 4 dims of each sample) ---"
    do j = 1, 8
      write(*,'(a,i0,a,4f8.4)') "  sample ", j, " : ", samples(1:4, j)
    end do
  endif

  !> SPX 交叉サンプリング
  call monolis_opt_vae_generate_spx(net, X, 8, samples)
  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,'(a)') " --- SPX-crossover samples (first 4 dims) ---"
    do j = 1, 8
      write(*,'(a,i0,a,4f8.4)') "  sample ", j, " : ", samples(1:4, j)
    end do
  endif

  call monolis_opt_vae_finalize(net)
  call monolis_global_finalize()
end program main
