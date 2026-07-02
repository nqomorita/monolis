!> GEMM カーネル最適化のマイクロベンチ
!> forward: pre(o,j) = sum_i W(i,o)*a_in(i,j)
!> 各方式を同一データで比較 (out_dim x nB 出力, in_dim 縮約)
program bench_gemm
  implicit none
  integer, parameter :: sp = 4
  integer, parameter :: in_dim = 256, out_dim = 512, nB = 256
  integer, parameter :: NITER = 2000
  real(sp), allocatable :: W(:,:), a_in(:,:), pre(:,:), ref(:,:)
  integer :: i, j, o, it, ts
  real(sp) :: s
  real(8) :: t0, t1
  real(sp) :: err

  allocate(W(in_dim,out_dim), a_in(in_dim,nB), pre(out_dim,nB), ref(out_dim,nB))
  call random_number(W); call random_number(a_in)
  W = W - 0.5_sp; a_in = a_in - 0.5_sp

  ! reference (CPU)
  do j = 1, nB
    do o = 1, out_dim
      s = 0.0_sp
      do i = 1, in_dim
        s = s + W(i,o)*a_in(i,j)
      end do
      ref(o,j) = s
    end do
  end do

  !$acc data copyin(W, a_in) create(pre)

  ! ---- v0: naive (compiler auto = gang collapse(2), vector reduction) ----
  !$acc parallel loop collapse(2) present(W,a_in,pre) private(s)
  do j = 1, nB
    do o = 1, out_dim
      s = 0.0_sp
      do i = 1, in_dim
        s = s + W(i,o)*a_in(i,j)
      end do
      pre(o,j) = s
    end do
  end do
  !$acc update self(pre)
  call check('v0 warmup', pre, ref, in_dim, out_dim, nB)
  call sync_time(t0)
  do it = 1, NITER
    !$acc parallel loop collapse(2) present(W,a_in,pre) private(s)
    do j = 1, nB
      do o = 1, out_dim
        s = 0.0_sp
        do i = 1, in_dim
          s = s + W(i,o)*a_in(i,j)
        end do
        pre(o,j) = s
      end do
    end do
  end do
  call sync_time(t1)
  write(*,'(a,f8.2,a)') 'v0 naive(auto)         : ', (t1-t0)*1.0d3, ' ms'

  ! ---- v1: tile(16,8) ----
  call sync_time(t0)
  do it = 1, NITER
    !$acc parallel loop tile(16,8) present(W,a_in,pre) private(s)
    do j = 1, nB
      do o = 1, out_dim
        s = 0.0_sp
        !$acc loop seq
        do i = 1, in_dim
          s = s + W(i,o)*a_in(i,j)
        end do
        pre(o,j) = s
      end do
    end do
  end do
  call sync_time(t1)
  !$acc update self(pre)
  call check('v1 tile16x8', pre, ref, in_dim, out_dim, nB)
  write(*,'(a,f8.2,a)') 'v1 tile(16,8)          : ', (t1-t0)*1.0d3, ' ms'

  ! ---- v2: tile(16,16) ----
  call sync_time(t0)
  do it = 1, NITER
    !$acc parallel loop tile(16,16) present(W,a_in,pre) private(s)
    do j = 1, nB
      do o = 1, out_dim
        s = 0.0_sp
        !$acc loop seq
        do i = 1, in_dim
          s = s + W(i,o)*a_in(i,j)
        end do
        pre(o,j) = s
      end do
    end do
  end do
  call sync_time(t1)
  !$acc update self(pre)
  call check('v2 tile16', pre, ref, in_dim, out_dim, nB)
  write(*,'(a,f8.2,a)') 'v2 tile(16,16) + seq   : ', (t1-t0)*1.0d3, ' ms'

  ! ---- v3: tile(32,16) ----
  call sync_time(t0)
  do it = 1, NITER
    !$acc parallel loop tile(32,16) present(W,a_in,pre) private(s)
    do j = 1, nB
      do o = 1, out_dim
        s = 0.0_sp
        !$acc loop seq
        do i = 1, in_dim
          s = s + W(i,o)*a_in(i,j)
        end do
        pre(o,j) = s
      end do
    end do
  end do
  call sync_time(t1)
  !$acc update self(pre)
  call check('v3 tile32x16', pre, ref, in_dim, out_dim, nB)
  write(*,'(a,f8.2,a)') 'v3 tile(32,16)         : ', (t1-t0)*1.0d3, ' ms'

  ! ---- v4: tile(16,32) ----
  call sync_time(t0)
  do it = 1, NITER
    !$acc parallel loop tile(16,32) present(W,a_in,pre) private(s)
    do j = 1, nB
      do o = 1, out_dim
        s = 0.0_sp
        !$acc loop seq
        do i = 1, in_dim
          s = s + W(i,o)*a_in(i,j)
        end do
        pre(o,j) = s
      end do
    end do
  end do
  call sync_time(t1)
  !$acc update self(pre)
  call check('v4 tile16x32', pre, ref, in_dim, out_dim, nB)
  write(*,'(a,f8.2,a)') 'v4 tile(16,32)         : ', (t1-t0)*1.0d3, ' ms'

  !$acc end data
contains
  subroutine sync_time(t)
    real(8), intent(out) :: t
    !$acc wait
    call cpu_time_wall(t)
  end subroutine
  subroutine cpu_time_wall(t)
    real(8), intent(out) :: t
    integer(8) :: c, r
    call system_clock(c, r)
    t = real(c,8)/real(r,8)
  end subroutine
  subroutine check(tag, p, rf, id, od, nb)
    character(*), intent(in) :: tag
    integer, intent(in) :: id, od, nb
    real(sp), intent(in) :: p(od,nb), rf(od,nb)
    real(sp) :: e
    e = maxval(abs(p - rf))
    write(*,'(a,a,es10.3)') tag, ' maxerr=', e
  end subroutine
end program bench_gemm
