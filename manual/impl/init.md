
## 初期化と終了処理

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

## monolis 構造体の初期化と終了処理

monolis を利用するためには `monolis_structure` 構造体を定義する。
`monolis_structure` 構造体は、行列の情報、求解パラメータ、並列計算のための通信情報を保持している。
`monolis_structure` 構造体は、これら情報を含んだ疎行列の構造体と理解するとよい。
`monolis_structure` 構造体をもつ変数は、`monolis_initialize` 関数による初期化の後に利用できる。

初期化処理は `monolis_initialize` 関数を用いて行う。
第1引数は `monolis_structure` 構造体である。
第2引数は入力ファイルが存在するディレクトリを指定する文字列である。
`monolis_global_initialize` 関数の実行後に利用できる。

終了処理は `monolis_finalize` を用いて行う。引数は `monolis_structure` 構造体である。

```fortran
subroutine sample
  use mod_monolis
  implicit none
  type(monolis_structure) :: A !> Ax = b の係数行列

  !> 初期化
  call monolis_global_initialize()
  call monolis_initialize(A, "./")

  !> 終了処理
  call monolis_finalize(A)
  call monolis_global_finalize()
end subroutine
```

