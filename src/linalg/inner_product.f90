!> ベクトル内積関数群
module mod_monolis_inner_product
  use mod_monolis_utils
  use mod_monolis_def_struc
  implicit none

contains

  !> @ingroup linalg
  !> ベクトル内積（整数型）
  subroutine monolis_inner_product_I(monolis, monoCOM, n_dof, X, Y, sum)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] ベクトル 1
    integer(kint), intent(in) :: X(:)
    !> [in] ベクトル 2
    integer(kint), intent(in) :: Y(:)
    !> [out] 内積結果
    integer(kint), intent(out) :: sum
    integer(kint) :: N
    real(kdouble) :: tdotp, tcomm

    N = monolis%MAT%N
    if(monoCOM%comm_size > 1) N = monoCOM%n_internal_vertex

    call monolis_inner_product_main_I(monoCOM, N*n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_I

  !> @ingroup linalg
  !> ベクトル内積（整数型、任意のベクトルサイズ）
  subroutine monolis_inner_product_V_I(monoCOM, m, X, Y, sum)
    implicit none
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in] 内部計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] ベクトル 1
    integer(kint), intent(in) :: X(:)
    !> [in] ベクトル 2
    integer(kint), intent(in) :: Y(:)
    !> [out] 内積結果
    integer(kint), intent(out) :: sum
    real(kdouble) :: tdotp, tcomm

    call monolis_inner_product_main_I(monoCOM, m, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_V_I

  !> @ingroup dev_linalg
  !> ベクトル内積（整数型、メイン関数）
  subroutine monolis_inner_product_main_I(monoCOM, m, X, Y, sum, tdotp, tcomm)
    implicit none
    !> [in] monoCOM 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 内部計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] ベクトル 1
    integer(kint), intent(in) :: X(:)
    !> [in] ベクトル 2
    integer(kint), intent(in) :: Y(:)
    !> [out] 内積結果
    integer(kint), intent(out) :: sum
    !> [inout] 内積時間
    real(kdouble), optional, intent(inout) :: tdotp
    !> [inout] 通信時間
    real(kdouble), optional, intent(inout) :: tcomm
    integer(kint) :: i
    real(kdouble) :: t1, t2, t3

    call monolis_std_debug_log_header("monolis_inner_product_main_I")

    t1 = monolis_get_time()
    sum = 0
!$omp parallel default(none) &
!$omp & shared(X, Y, sum) &
!$omp & firstprivate(m) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, m
      sum = sum + X(i)*Y(i)
    enddo
!$omp end do
!$omp end parallel

    t2 = monolis_get_time()
    call monolis_allreduce_I1(sum, monolis_mpi_sum, monoCOM%comm)
    t3 = monolis_get_time()

    if(present(tdotp)) tdotp = tdotp + t3 - t1
    if(present(tcomm)) tcomm = tcomm + t3 - t2
  end subroutine monolis_inner_product_main_I

  !> @ingroup linalg
  !> ベクトル内積（実数型）
  subroutine monolis_inner_product_R(monolis, monoCOM, n_dof, X, Y, sum)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] ベクトル 1
    real(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    real(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    real(kdouble), intent(out) :: sum
    integer(kint) :: N
    real(kdouble) :: tdotp, tcomm

    N = monolis%MAT%N
    if(monoCOM%comm_size > 1) N = monoCOM%n_internal_vertex

    call monolis_inner_product_main_R(monoCOM, N*n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_R

  !> @ingroup linalg
  !> ベクトル内積（実数型、任意のベクトルサイズ）
  subroutine monolis_inner_product_V_R(monoCOM, m, X, Y, sum)
    implicit none
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in] 内部計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] ベクトル 1
    real(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    real(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    real(kdouble), intent(out) :: sum
    real(kdouble) :: tdotp, tcomm

    call monolis_inner_product_main_R(monoCOM, m, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_V_R

  !> @ingroup dev_linalg
  !> ベクトル内積（実数型、メイン関数）
  subroutine monolis_inner_product_main_R(monoCOM, m, X, Y, sum, tdotp, tcomm)
    implicit none
    !> [in] monoCOM 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 内部計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] ベクトル 1
    real(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    real(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    real(kdouble), intent(out) :: sum
    !> [inout] 内積時間
    real(kdouble), optional, intent(inout) :: tdotp
    !> [inout] 通信時間
    real(kdouble), optional, intent(inout) :: tcomm
    integer(kint) :: i
    real(kdouble) :: t1, t2, t3

    call monolis_std_debug_log_header("monolis_inner_product_main_R")

    t1 = monolis_get_time()
    sum = 0.0d0
!$omp parallel default(none) &
!$omp & shared(X, Y, sum) &
!$omp & firstprivate(m) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, m
      sum = sum + X(i)*Y(i)
    enddo
!$omp end do
!$omp end parallel

    t2 = monolis_get_time()
    call monolis_allreduce_R1(sum, monolis_mpi_sum, monoCOM%comm)
    t3 = monolis_get_time()

    if(present(tdotp)) tdotp = tdotp + t3 - t1
    if(present(tcomm)) tcomm = tcomm + t3 - t2
  end subroutine monolis_inner_product_main_R

  !> @ingroup linalg
  !> ベクトル内積（擬似四倍精度実数型）
  subroutine monolis_inner_product_R_N128(monolis, monoCOM, n_dof, X, Y, sum)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] ベクトル 1
    real(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    real(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    real(kdouble), intent(out) :: sum
    integer(kint) :: N
    real(kdouble) :: tdotp, tcomm

    N = monolis%MAT%N
    if(monoCOM%comm_size > 1) N = monoCOM%n_internal_vertex

    call monolis_inner_product_main_R_N128(monoCOM, N, n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_R_N128

  !> @ingroup dev_linalg
  !> ベクトル内積（擬似四倍精度実数型、メイン関数）
  subroutine monolis_inner_product_main_R_N128(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
    implicit none
    !> [in] monoCOM 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 内部計算点数
    integer(kint), intent(in) :: n
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] ベクトル 1
    real(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    real(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    real(kdouble), intent(out) :: sum
    !> [inout] 内積時間
    real(kdouble), optional, intent(inout) :: tdotp
    !> [inout] 通信時間
    real(kdouble), optional, intent(inout) :: tcomm
    type(monolis_R_N128) :: sum_N128, a, b
    integer(kint) :: i
    real(kdouble) :: t1, t2, t3

    call monolis_std_debug_log_header("monolis_inner_product_main_R_N128")

    t1 = monolis_get_time()
    sum_N128 = monolis_conv_R_to_R_N128(0.0d0)

    do i = 1, n * n_dof
      a = monolis_conv_R_to_R_N128(X(i)*Y(i))
      b = monolis_copy_R_N128(sum_N128)
      call monolis_add_R_N128(a, b, sum_N128)
    enddo

    t2 = monolis_get_time()
    call monolis_allreduce_R1_N128(sum_N128, monolis_mpi_sum, monoCOM%comm)
    t3 = monolis_get_time()

    sum = monolis_conv_R_N128_to_R(sum_N128)

    if(present(tdotp)) tdotp = tdotp + t3 - t1
    if(present(tcomm)) tcomm = tcomm + t3 - t2
  end subroutine monolis_inner_product_main_R_N128

  !> @ingroup linalg
  !> ベクトル内積（擬似四倍精度実数型、ランク 0 に全ての配列要素を集め、絶対値昇順で内積計算）
  subroutine monolis_global_sorted_inner_product_main_R_N128(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
    implicit none
    !> [in] monoCOM 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 内部計算点数
    integer(kint), intent(in) :: n
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] ベクトル 1
    real(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    real(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    real(kdouble), intent(out) :: sum
    !> [inout] 内積時間
    real(kdouble), intent(inout) :: tdotp
    !> [inout] 通信時間
    real(kdouble), intent(inout) :: tcomm
    integer(kint) :: n_global, i, my_rank, comm_size, root
    real(kdouble) :: sum_global
    integer(kint), allocatable :: rc(:), disp(:)
    real(kdouble), allocatable :: X_global(:)
    real(kdouble), allocatable :: Y_global(:)
    real(kdouble), allocatable :: XY_global(:)
    real(kdouble), allocatable :: XY_abs_global(:)
    type(monolis_R_N128) :: sum_N128, a, b

    !> gather part
    my_rank = monolis_mpi_get_local_my_rank(monoCOM%comm)
    comm_size = monolis_mpi_get_local_comm_size(monoCOM%comm)

    n_global = n
    call monolis_allreduce_I1(n_global, monolis_mpi_sum, monoCOM%comm)

    !if(my_rank == 0)then
      call monolis_alloc_R_1d(X_global, n_global*n_dof)
      call monolis_alloc_R_1d(Y_global, n_global*n_dof)
      call monolis_alloc_R_1d(XY_global, n_global*n_dof)
      call monolis_alloc_R_1d(XY_abs_global, n_global*n_dof)
    !endif

    call monolis_alloc_I_1d(rc, comm_size)
    call monolis_alloc_I_1d(disp, comm_size)

    call monolis_allgather_I1(n*n_dof, rc, monoCOM%comm)
    do i = 1, comm_size - 1
      disp(i + 1) = disp(i) + rc(i)
    enddo

    root = 0
    call monolis_gather_V_R(X, n*n_dof, X_global, rc, disp, root, monoCOM%comm)
    call monolis_gather_V_R(Y, n*n_dof, Y_global, rc, disp, root, monoCOM%comm)

    !> main inner product routine
    sum_global = 0.0d0

    if(my_rank == 0)then
      do i = 1, n_global*n_dof
        XY_global(i) = X_global(i)*Y_global(i)
      enddo

      XY_abs_global = dabs(XY_global)
      call monolis_qsort_R_2d(XY_abs_global, XY_global, 1, n_global*n_dof)

      sum_N128 = monolis_conv_R_to_R_N128(0.0d0)
      do i = 1, n_global*n_dof
        a = monolis_conv_R_to_R_N128(XY_global(i))
        b = monolis_copy_R_N128(sum_N128)
        call monolis_add_R_N128(a, b, sum_N128)
      enddo
      sum_global = monolis_conv_R_N128_to_R(sum_N128)
    endif

    !sum = sum_global
    call monolis_allreduce_R1(sum_global, monolis_mpi_sum, monoCOM%comm)
    sum = sum_global
  end subroutine monolis_global_sorted_inner_product_main_R_N128

  !> @ingroup linalg
  !> ベクトル内積（複素数型）
  subroutine monolis_inner_product_C(monolis, monoCOM, n_dof, X, Y, sum)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] ベクトル 1
    complex(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    complex(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    complex(kdouble), intent(out) :: sum
    integer(kint) :: N
    real(kdouble) :: tdotp, tcomm

    N = monolis%MAT%N
    if(monoCOM%comm_size > 1) N = monoCOM%n_internal_vertex

    call monolis_inner_product_main_C(monoCOM, N*n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_C

  !> @ingroup linalg
  !> ベクトル内積（複素数型、任意のベクトルサイズ）
  subroutine monolis_inner_product_V_C(monoCOM, m, X, Y, sum)
    implicit none
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in] 内部計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] ベクトル 1
    complex(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    complex(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    complex(kdouble), intent(out) :: sum
    real(kdouble) :: tdotp, tcomm

    call monolis_inner_product_main_C(monoCOM, m, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_V_C

  !> @ingroup dev_linalg
  !> ベクトル内積（複素数型、メイン関数）
  subroutine monolis_inner_product_main_C(monoCOM, m, X, Y, sum, tdotp, tcomm)
    implicit none
    !> [in] monoCOM 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 内部計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] ベクトル 1
    complex(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    complex(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    complex(kdouble), intent(out) :: sum
    !> [inout] 内積時間
    real(kdouble), optional, intent(inout) :: tdotp
    !> [inout] 通信時間
    real(kdouble), optional, intent(inout) :: tcomm
    integer(kint) :: i
    real(kdouble) :: t1, t2, t3

    call monolis_std_debug_log_header("monolis_inner_product_main_C")

    t1 = monolis_get_time()
    sum = 0.0d0
!$omp parallel default(none) &
!$omp & shared(X, Y, sum) &
!$omp & firstprivate(m) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, m
      sum = sum + X(i)*Y(i)
    enddo
!$omp end do
!$omp end parallel

    t2 = monolis_get_time()
    call monolis_allreduce_C1(sum, monolis_mpi_sum, monoCOM%comm)
    t3 = monolis_get_time()

    if(present(tdotp)) tdotp = tdotp + t3 - t1
    if(present(tcomm)) tcomm = tcomm + t3 - t2
  end subroutine monolis_inner_product_main_C

  !> @ingroup dev_linalg
  !> ベクトル内積（実数型、メイン関数、通信なし）
  subroutine monolis_inner_product_main_R_no_comm(m, X, Y, sum)
    implicit none
    !> [in] 内部計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] ベクトル 1
    real(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    real(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    real(kdouble), intent(out) :: sum
    integer(kint) :: i

    sum = 0.0d0
!$omp parallel default(none) &
!$omp & shared(X, Y, sum) &
!$omp & firstprivate(m) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, m
      sum = sum + X(i)*Y(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_inner_product_main_R_no_comm
end module mod_monolis_inner_product