
### 疎行列の非零パターンの決定

#### 要素情報からの決定（要素種類が単一の場合）

この関数は要素から非零パターンの決定する関数である。
要素種類が単一の場合の疎行列の非零パターンの取得は、`monolis_get_nonzero_pattern` 関数で行う。
`monolis_global_initialize` 関数および `monolis_initialize` 関数の実行後に利用できる。
以下に図示する 2 節点トラス 3 要素からなる行列における例を示す。

<img src="./fig/nonzero.svg" height=300px>

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
  call monolis_initialize(A, "./")

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

#### 要素情報からの決定（要素種類が複数の場合）

この関数は要素から非零パターンの決定する関数である。
要素種類が複数の場合の疎行列の非零パターンの取得は、`monolis_get_nonzero_pattern_by_connectivity` 関数で行う。
`monolis_global_initialize` 関数および `monolis_initialize` 関数の実行後に利用できる。
以下に図示する 2 節点トラス 3 要素からなる行列における例を示す。

<img src="./fig/nonzero.svg" height=300px>

```fortran
program main
  use mod_monolis
  implicit none
  type(monolis_structure) :: A !> Ax = b の係数行列
  integer(kint) :: nnode !> 節点数
  integer(kint) :: nelem !> 要素数
  integer(kint) :: index(:) !> 要素定義の index 配列
  integer(kint) :: item(:) !> 要素定義の item 配列
  integer(kint) :: ndof !> 1 節点あたりの自由度

  !> 初期化
  call monolis_global_initialize()
  call monolis_initialize(A, "./")

  !> 疎行列の非零パターンの決定
  nnode = 4
  nelem = 3
  ndof = 2

  index(0) = 0
  index(1) = 2
  index(2) = 4
  index(3) = 6
  index(4) = 8
  item(1) = 1; item(2) = 3 !> 要素番号 1 の定義
  item(3) = 3; item(4) = 2 !> 要素番号 2 の定義
  item(5) = 2; item(6) = 4 !> 要素番号 3 の定義

  call monolis_get_nonzero_pattern_by_connectivity(A, nnode, ndof, nelem, index, item)

  !> 終了処理
  call monolis_finalize(A)
  call monolis_global_finalize()
end program main
```

#### 節点グラフからの決定

この関数は節点グラフから非零パターンの決定する関数である。
疎行列の非零パターンの取得は、`monolis_get_nonzero_pattern_by_nodal` 関数で行う。
`monolis_global_initialize` 関数および `monolis_initialize` 関数の実行後に利用できる。
以下に図示する 2 節点トラス 3 要素からなる行列における例を示す。

<img src="./fig/nonzero.svg" height=300px>

```fortran
program main
  use mod_monolis
  implicit none
  type(monolis_structure) :: A !> Ax = b の係数行列
  integer(kint) :: nnode !> 節点数
  integer(kint) :: index(:) !> 要素定義の index 配列
  integer(kint) :: item(:) !> 要素定義の item 配列
  integer(kint) :: ndof !> 1 節点あたりの自由度

  !> 初期化
  call monolis_global_initialize()
  call monolis_initialize(A, "./")

  !> 疎行列の非零パターンの決定
  nnode = 4
  ndof = 2

  index(0) = 0
  index(1) = 1
  index(2) = 3
  index(3) = 5
  index(4) = 6
  !> 自身の節点番号は含まない
  item(1) = 3 !> 節点 1 と接続している節点番号
  item(2) = 3 !> 節点 2 と接続している節点番号
  item(3) = 4 !> 節点 2 と接続している節点番号
  item(4) = 1 !> 節点 3 と接続している節点番号
  item(5) = 2 !> 節点 3 と接続している節点番号
  item(6) = 2 !> 節点 4 と接続している節点番号
  call monolis_get_nonzero_pattern_by_nodal(A, nnode, ndof, index, item)

  !> 終了処理
  call monolis_finalize(A)
  call monolis_global_finalize()
end program main
```


### 行列値の初期化

疎行列パターンを保持したまま行列の値を 0 初期化する場合は、`monolis_clear_mat_value` 関数で行う。
この関数では、行列、右辺ベクトル、解ベクトルを全て 0 に初期化する。

```fortran
  type(monolis_structure) :: A !> y = Ax の係数行列

  !> 行列値の初期化
  call monolis_clear_mat_value(A)
```



