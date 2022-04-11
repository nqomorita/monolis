

## 全体剛性行列への足し込み

### 要素剛性行列の足し込み

要素剛性行列を全体剛性行列へ足し込む場合は、`monolis_add_matrix_to_sparse_matrix` 関数で行う。
`monolis_get_nonzero_pattern` 関数の実行後に利用できる。

```fortran
  type(monolis_structure) :: A !> Ax = b の係数行列
  integer(kint) :: nbase_func !> 要素の次数（要素を構成する節点数）
  integer(kint) :: connectivity(2) !> 要素を構成する節点番号
  real(kdouble) :: local_mat(4,4) !> 要素剛性行列
                                  !> local_mat(nbase_func*ndof, nbase_func*ndof)

  !> 疎行列への足し込み
  nbase_func = 2
  connectivity(1) = 1; connectivity(2) = 3 !> 要素番号 1 の定義
  call get_localmat(local_mat) !> user subroutine
  call monolis_add_matrix_to_sparse_matrix(A, nbase_func, connectivity, local_mat)
```

### 要素剛性行列の足し込み（副対角部分）

要素剛性行列を全体剛性行列の副対角部分へ足し込む場合は、`monolis_add_matrix_to_sparse_matrix_offdiag` 関数で行う。
`monolis_get_nonzero_pattern` 関数の実行後に利用できる。

```fortran
  type(monolis_structure) :: A !> Ax = b の係数行列
  integer(kint) :: n1 !> 行列成分を足し込む行番号の個数
  integer(kint) :: n2 !> 行列成分を足し込む列番号の個数
  integer(kint) :: c1(2) !> 要素を構成する行番号
  integer(kint) :: c2(2) !> 要素を構成する列番号
  real(kdouble) :: local_mat(4,4) !> 要素剛性行列
                                  !> local_mat(n1*ndof, n2*ndof)

  !> 疎行列への足し込み
  n1 = 2
  n2 = 2
  c1(1) = 1; c1(2) = 3 !> 要素を構成する行番号
  c2(1) = 3; c2(2) = 4 !> 要素を構成する列番号
  call get_localmat(local_mat) !> user subroutine
  call monolis_add_matrix_to_sparse_matrix(A, nbase_func, n1, n2, c1, c2, local_mat)
```

### 行列成分の足し込み

行列成分を全体剛性行列へ足し込む場合は、`monolis_add_scalar_to_sparse_matrix` 関数で行う。
`monolis_get_nonzero_pattern` 関数の実行後に利用できる。

```fortran
  type(monolis_structure) :: A !> Ax = b の係数行列
  integer(kint) :: i !> 行列成分を足し込む行番号
  integer(kint) :: j !> 行列成分を足し込む列番号
  integer(kint) :: sub_i !> 行列成分を足し込む小行列の行番号
  integer(kint) :: sub_j !> 行列成分を足し込む小行列の列番号
  real(kdouble) :: val !> 要素剛性行列
  !> (i, j) 番目に位置するブロック行列の (sub_i, sub_j) 自由度に値 val を足し込む

  !> 疎行列への足し込み
  do !> 行列成分の個数
    call monolis_add_scalar_to_sparse_matrix(A, i, j, sub_i, sub_j, val)
  enddo
```

