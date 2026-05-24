!> BLAS インターフェース（fact 内部用）
!> @note monolis 全体は openblas / Intel MKL / Apple Accelerate などに直接リンク済み。
!>       fact では BLAS の標準名を external 宣言してそのまま呼ぶ。
module mod_monolis_fact_blas
  use mod_monolis_utils
  implicit none

  interface
    subroutine daxpy(n, da, dx, incx, dy, incy)
      import :: kint, kdouble
      integer(kint), intent(in) :: n, incx, incy
      real(kdouble), intent(in) :: da
      real(kdouble), intent(in) :: dx(*)
      real(kdouble), intent(inout) :: dy(*)
    end subroutine

    subroutine dscal(n, da, dx, incx)
      import :: kint, kdouble
      integer(kint), intent(in) :: n, incx
      real(kdouble), intent(in) :: da
      real(kdouble), intent(inout) :: dx(*)
    end subroutine

    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: kint, kdouble
      character(len=1), intent(in) :: trans
      integer(kint), intent(in) :: m, n, lda, incx, incy
      real(kdouble), intent(in) :: alpha, beta
      real(kdouble), intent(in) :: a(lda, *), x(*)
      real(kdouble), intent(inout) :: y(*)
    end subroutine

    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: kint, kdouble
      character(len=1), intent(in) :: transa, transb
      integer(kint), intent(in) :: m, n, k, lda, ldb, ldc
      real(kdouble), intent(in) :: alpha, beta
      real(kdouble), intent(in) :: a(lda, *), b(ldb, *)
      real(kdouble), intent(inout) :: c(ldc, *)
    end subroutine

    subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: kint, kdouble
      character(len=1), intent(in) :: side, uplo, transa, diag
      integer(kint), intent(in) :: m, n, lda, ldb
      real(kdouble), intent(in) :: alpha
      real(kdouble), intent(in) :: a(lda, *)
      real(kdouble), intent(inout) :: b(ldb, *)
    end subroutine

    subroutine dlaset(uplo, m, n, alpha, beta, a, lda)
      import :: kint, kdouble
      character(len=1), intent(in) :: uplo
      integer(kint), intent(in) :: m, n, lda
      real(kdouble), intent(in) :: alpha, beta
      real(kdouble), intent(inout) :: a(lda, *)
    end subroutine
  end interface
end module mod_monolis_fact_blas
