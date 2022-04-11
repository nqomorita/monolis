
### 行列ベクトル積

行列ベクトル積の計算は、`monolis_matvec_product` 関数で行う。

```fortran
  type(monolis_structure) :: A !> y = Ax の係数行列
  real(kdouble), allocatable :: x(:) !> 右辺ベクトル
  real(kdouble), allocatable :: y(:) !> 解ベクトル

  !> 行列ベクトル積
  allocate(x(nnode*ndof), source = 0.0d0)
  allocate(y(nnode*ndof), source = 0.0d0)
  call get_RHS(x) !> user subroutine
  call monolis_matvec_product(A, x, y)
```



### ベクトル内積

ベクトル内積の計算は、`monolis_inner_product` 関数で行う。


