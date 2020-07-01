# monolis

- Morita's non-overlapping / overlapping DDM based linear equation solver

## 既存ソルバへの組み込み

### 初期化と終了処理

全体初期化処理は `monolis_global_initialize` 関数を用いて行う。引数なし。

全体終了処理は `monolis_global_finalize` を用いて行う。引数なし。

全体初期化処理と全体終了処理は、プログラム実行中にそれぞれ 1 回しか実行できない。

```fortran
program main
  use mod_monolis
  implicit none

  !> 初期化
  call monolis_global_initialize()

  !> 終了処理
  call monolis_global_finalize()
end program main
```

### 時間計測

時間計測は `monolis_get_time` 関数を用いて行う。引数なし。戻り値は倍精度浮動小数点数型である。

```fortran
program main
  use mod_monolis
  implicit none
  real(8) :: t1, t2

  call monolis_global_initialize()

  !> 時間計測
  t1 = monolis_get_time()

  !> 時間計測
  t2 = monolis_get_time()

  !> 計測時間の出力
  write(*,*)"** monolis time", t2 - t1

  call monolis_global_finalize()
end program main
```

### monolis 構造体

monolis を利用するためには `monolis_structure` 構造体を定義する。
`monolis_structure` 構造体は、行列の情報、求解パラメータ、並列計算のための通信情報を保持している。
`monolis_structure` 構造体は、これら情報を含んだ疎行列の構造体と理解するとよい。
`monolis_structure` 構造体をもつ変数は、`monolis_initialize` 関数による初期化の後に利用できる。

初期化処理は `monolis_initialize` 関数を用いて行う。引数は `monolis_structure` 構造体である。
`monolis_global_initialize` 関数の実行後に利用できる。

終了処理は `monolis_finalize` を用いて行う。引数は `monolis_structure` 構造体である。

```fortran
subroutine sample
  use mod_monolis
  implicit none
  type(monolis_structure) :: A !> Ax = b の係数行列

  !> 初期化
  call monolis_global_initialize()
  call monolis_initialize(A)

  !> 終了処理
  call monolis_finalize(A)
  call monolis_global_finalize()
end subroutine
```

### 疎行列の非零パターンの決定

疎行列の非零パターンの取得は、`monolis_get_nonzero_pattern` 関数で行う。
`monolis_global_initialize` 関数および `monolis_initialize` 関数の実行後に利用できる。
以下に図示する 2 節点トラス 3 要素からなる行列における例を示す。

<img src="./nonzero.svg" height=300px>

```fortran
program main
  use mod_monolis
  implicit none
  type(monolis_structure) :: A !> Ax = b の係数行列
  integer(kint) :: nnode !> 節点数
  integer(kint) :: nelem !> 要素数
  integer(kint) :: nbase_func !> 要素の次数（要素を構成する節点数）
  integer(kint) :: ndof !> 1 節点あたりの自由度
  integer(kint) :: elem(2,4) !> elem(nbase_func, nnode) 要素定義（例：2 節点トラス、3 要素）

  !> 初期化
  call monolis_global_initialize()
  call monolis_initialize(A)

  !> 疎行列の非零パターンの決定
  nnode = 4
  nelem = 3
  nbase_func = 2
  ndof = 2
  elem(1,1) = 1; elem(2,1) = 3 !> 要素番号 1 の定義
  elem(1,2) = 3; elem(2,2) = 2 !> 要素番号 2 の定義
  elem(1,3) = 2; elem(2,3) = 4 !> 要素番号 3 の定義

  call monolis_get_nonzero_pattern(A, nnode, nbase_func, ndof, nelem, elem)

  !> 終了処理
  call monolis_finalize(A)
  call monolis_global_finalize()
end program main
```

#### 要素行列への要素行列の足し込み

要素行列への行列成分の足し込みは、`monolis_assemble_sparse_matrix` 関数で行う。

```fortran
  type(monolis_structure) :: A !> Ax = b の係数行列
  integer(kint) :: nbase_func !> 要素の次数（要素を構成する節点数）
  integer(kint) :: connectivity(2) !> 要素を構成する節点番号
  real(kdouble) :: local_mat(4,4) !> local_mat(nbase_func*ndof, nbase_func*ndof)、要素剛性行列

  !> 疎行列への足し込み
  nbase_func = 2
  connectivity(1) = 1; connectivity(2) = 3 !> 要素番号 1 の定義
  call get_localmat(local_mat) !> for example
  call monolis_assemble_sparse_matrix(A, nbase_func, connectivity, local_mat)
```

#### Dirichlet 境界条件の追加

要素行列への Dirichlet 境界条件の追加は、`monolis_set_Dirichlet_bc` 関数で行う。

```fortran
  type(monolis_structure) :: A !> Ax = b の係数行列
  integer(kint) :: node_id !> 節点番号
  integer(kint) :: ndof_bc !> 境界条件を与える自由度番号（[1, ndof] の値をとる）
  real(kdouble) :: val !> 境界条件の設定値

  !> Dirichlet 境界条件の追加
  node_id = 3
  ndof_bc = 2
  val = 0.0d0

  !> 節点番号 3 の第 2 自由度成分の値を 0.0 に固定
  call monolis_set_Dirichlet_bc(A, node_id, ndof_bc, val)
```

<!--
#### Neumann 境界条件の設定

右辺ベクトルへのベクトル成分の足し込みは、`monolis_assemble_RHS` 関数で行う。

```fortran
  type(monolis_structure) :: A !> Ax = b の係数行列
  integer(kint) :: node_id !> 節点番号
  integer(kint) :: ndof_bc !> 境界条件を与える自由度番号（[1, ndof] の値をとる）
  real(kdouble) :: val !> 境界条件の設定値

  !> Neumann 境界条件の追加
  node_id = 1
  ndof_bc = 2
  val = 1.0d0

  !> 節点番号 1 の第 2 自由度成分に対応する右辺値に 1.0 を加算
  call monolis_assemble_RHS(A, node_id, ndof_bc, val)
```
-->

#### 方程式の求解

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
  call get_RHS(b) !> for example
  call monolis_solve(A, b, x)
```

<!--
#### 解ベクトルの取得

解ベクトルの取得は、`monolis_get_solution` 関数で行う。

```fortran
  type(monolis_structure) :: A !> Ax = b の係数行列
  real(8), allocatable :: x(:)

  !> 解ベクトルの取得
  call monolis_get_solution(A, x)
```
-->




