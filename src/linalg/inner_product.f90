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

    call monolis_inner_product_main_I(monoCOM, N, n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_I

  !> @ingroup linalg
  !> ベクトル内積（整数型、任意のベクトルサイズ）
  subroutine monolis_inner_productV_I(monolis, monoCOM, n, n_dof, X, Y, sum)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in] 内部計算点数
    integer(kint), intent(in) :: n
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] ベクトル 1
    integer(kint), intent(in) :: X(:)
    !> [in] ベクトル 2
    integer(kint), intent(in) :: Y(:)
    !> [out] 内積結果
    integer(kint), intent(out) :: sum
    real(kdouble) :: tdotp, tcomm

    call monolis_inner_product_main_I(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_productV_I

  !> @ingroup dev_linalg
  !> ベクトル内積（整数型、メイン関数）
  subroutine monolis_inner_product_main_I(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
    implicit none
    !> [in] monoCOM 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 内部計算点数
    integer(kint), intent(in) :: n
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
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
!$omp & firstprivate(n, n_dof) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, n * n_dof
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

    call monolis_inner_product_main_R(monoCOM, N, n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_R

  !> @ingroup linalg
  !> ベクトル内積（実数型、任意のベクトルサイズ）
  subroutine monolis_inner_productV_R(monolis, monoCOM, n, n_dof, X, Y, sum)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
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
    real(kdouble) :: tdotp, tcomm

    call monolis_inner_product_main_R(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_productV_R

  !> @ingroup dev_linalg
  !> ベクトル内積（実数型、メイン関数）
  subroutine monolis_inner_product_main_R(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
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
    integer(kint) :: i
    real(kdouble) :: t1, t2, t3

    call monolis_std_debug_log_header("monolis_inner_product_main_R")

    t1 = monolis_get_time()
    sum = 0.0d0
!$omp parallel default(none) &
!$omp & shared(X, Y, sum) &
!$omp & firstprivate(n, n_dof) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, n * n_dof
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

    call monolis_inner_product_main_C(monoCOM, N, n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_C

  !> @ingroup linalg
  !> ベクトル内積（複素数型、任意のベクトルサイズ）
  subroutine monolis_inner_productV_C(monolis, monoCOM, n, n_dof, X, Y, sum)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in] 内部計算点数
    integer(kint), intent(in) :: n
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] ベクトル 1
    complex(kdouble), intent(in) :: X(:)
    !> [in] ベクトル 2
    complex(kdouble), intent(in) :: Y(:)
    !> [out] 内積結果
    complex(kdouble), intent(out) :: sum
    real(kdouble) :: tdotp, tcomm

    call monolis_inner_product_main_C(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_productV_C

  !> @ingroup dev_linalg
  !> ベクトル内積（複素数型、メイン関数）
  subroutine monolis_inner_product_main_C(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
    implicit none
    !> [in] monoCOM 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 内部計算点数
    integer(kint), intent(in) :: n
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
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
!$omp & firstprivate(n, n_dof) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, n * n_dof
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
  subroutine monolis_inner_product_main_R_no_comm(n, n_dof, X, Y, sum)
    implicit none
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
    integer(kint) :: i

    sum = 0.0d0
!$omp parallel default(none) &
!$omp & shared(X, Y, sum) &
!$omp & firstprivate(n, n_dof) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, n * n_dof
      sum = sum + X(i)*Y(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_inner_product_main_R_no_comm
end module mod_monolis_inner_product