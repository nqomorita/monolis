
### 時間計測

時間計測は `monolis_get_time` 関数または `monolis_get_time_sync` 関数を用いて行う。引数なし。戻り値は倍精度浮動小数点数型である。
`monolis_get_time_sync` 関数は、並列計算実行時に全プロセスの同期処理（MPI_barrier）が実行される。

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



### 同期処理

並列計算時の同期処理は、`monolis_barrier` 関数で行う。

```fortran
  type(monolis_structure) :: A !> y = Ax の係数行列

  !> 並列計算の同期処理
  call monolis_barrier(A)
```


