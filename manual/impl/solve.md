

## Dirichlet 境界条件の追加

要素行列への Dirichlet 境界条件の追加は、`monolis_set_Dirichlet_bc` 関数で行う。
全体剛性行列を作成した後に利用できる。

```fortran
  type(monolis_structure) :: A !> Ax = b の係数行列
  integer(kint) :: node_id !> 節点番号
  integer(kint) :: ndof_bc !> 境界条件を与える自由度番号（[1, ndof] の値をとる）
  real(kdouble) :: b(:) !> 右辺ベクトル
  real(kdouble) :: val !> 境界条件の設定値

  !> Dirichlet 境界条件の追加
  node_id = 3
  ndof_bc = 2
  val = 0.0d0

  !> 節点番号 3 の第 2 自由度成分の値を 0.0 に固定
  call monolis_set_Dirichlet_bc(A, b, node_id, ndof_bc, val)
```


## 方程式の求解

解ベクトルの取得は、`monolis_solve` 関数で行う。
右辺ベクトルはユーザが独自で作成し、monolis に引数として渡す。

```fortran
  type(monolis_structure) :: A !> Ax = b の係数行列
  real(kdouble), allocatable :: x(:) !> 解ベクトル
  real(kdouble), allocatable :: b(:) !> 右辺ベクトル
  integer(kint) :: nnode !> 節点数
  integer(kint) :: ndof !> 1 節点あたりの自由度

  !> 解ベクトルの取得
  allocate(x(nnode*ndof), source = 0.0d0)
  allocate(b(nnode*ndof), source = 0.0d0)
  call get_RHS(b) !> user subroutine
  call monolis_solve(A, b, x)
```


