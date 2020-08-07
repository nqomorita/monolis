# monolis

- Morita's non-overlapping / overlapping DDM based linear equation solver

## 求解時のパラメータ設定

求解時のパラメータ設定は、`monolis_solve` 関数が呼ばれる直前までに行う。

### 関数一覧

- 反復解法の設定
    - monolis_set_method(monolis, method_param)
- 反復法前処理の設定
    - monolis_set_precond(monolis, prec_param)
- 最大反復回数の設定
    - monolis_set_maxiter(monolis, int_param)
- 収束判定閾値の設定
    - monolis_set_tolerance(monolis, double_param)

<!--
- 入力係数行列のスケーリングの設定
    - monolis_param_set_is_scaling(monolis, bool)
- 入力係数行列のリオーダリングの設定
    - monolis_param_set_is_reordering(monolis, bool)
- 入力解ベクトルの初期化の設定
    - monolis_param_set_is_init_x(monolis, bool)
- 入力係数行列を対称行列と見なすかの設定
    - monolis_param_set_is_sym_matrix(monolis, bool)
- デバッグログ出力の設定
    - monolis_param_set_is_debug(monolis, bool)
- 入力係数行列の対角成分に零成分が含まれの設定
    - monolis_param_set_is_check_diag(monolis, bool)
-
    - monolis_param_set_show_iterlog(monolis, bool)
    - monolis_param_set_show_time(monolis, bool)
    - monolis_param_set_show_summary(monolis, bool)
-->

### 反復解法の設定

```fortran
  monolis_set_method(
      monolis_structure   monolis,
      int                 method_param)
```

- デフォルト値：monolis_iter_CG
- method_param は以下のパラメータを設定する

| 変数名 | 説明 |
| ---- | ---- |
| monolis_iter_CG | Conjugate Gradient 法 |
| monolis_iter_GropCG | Grop's CG 法 |
| monolis_iter_PipeCG | Pipelined CG 法 |
| monolis_iter_PipeCR | Pipelined Conjugate Residual 法 |
| monolis_iter_BiCGSTAB | Bi-Conjugate Gradient Stabilization 法 |
| monolis_iter_PipeBiCGSTAB | Pipelined BiCGSTAB 法 |
| monolis_iter_BiCGSTAB_noprec | BiCGSTAB 法 (前処理無) |
| monolis_iter_CABiCGSTAB_noprec | Comunnication-avoiding BiCGSTAB 法 (前処理無) |
| monolis_iter_PipeBiCGSTAB_noprec | Pipelined BiCGSTAB 法 (前処理無) |
| monolis_iter_SOR | Successive Over-Relaxation 法 |
| monolis_iter_IR | Iterative Refinment 法 |

利用例

```fortran
  use mod_monolis
  implicit none
  type(monolis_structure) :: A !> Ax = b の係数行列
  call monolis_set_method(A, monolis_iter_CG) !> CG 法の設定
```

### 反復法前処理の設定

```fortran
  monolis_set_precond(
      monolis_structure   monolis,
      int                 prec_param)
```

- デフォルト値：monolis_prec_DIAG
- prec_param は以下のパラメータを設定する

| 変数名 | 説明 |
| ---- | ---- |
| monolis_prec_NONE | 前処理無 |
| monolis_prec_DIAG | 対角スケーリング前処理 |
| monolis_prec_ILU | Incomplete LU 分解前処理 |
| monolis_prec_JACOBI | Jacobi 前処理 |
| monolis_prec_SOR | SOR 前処理 |
| monolis_prec_SAINV | Stabilized Approximate Inverse 前処理 |
| monolis_prec_RIF | Robust Incomplete Factorization 前処理 |

<!--
| monolis_prec_SPIKE | SPIKE 前処理 |
| monolis_prec_DIRECT | 直接法 |
| monolis_prec_MUMPS | MUMPS |
-->

利用例

```fortran
  use mod_monolis
  implicit none
  type(monolis_structure) :: A !> Ax = b の係数行列
  call monolis_set_precond(A, monolis_prec_DIAG) !> 対角スケーリング前処理の設定
```

### 最大反復回数の設定

```fortran
  monolis_set_maxiter(
      monolis_structure   monolis,
      int                 maxiter)
```

- デフォルト値：1000

利用例

```fortran
  use mod_monolis
  implicit none
  type(monolis_structure) :: A !> Ax = b の係数行列
  call monolis_set_maxiter(A, 10000) !> 反復法の最大反復回数を 10000 回に設定
```

### 収束判定閾値の設定

```fortran
  monolis_set_tolerance(
      monolis_structure   monolis,
      double              threshold)
```

- デフォルト値：1.0e-8

利用例

```fortran
  use mod_monolis
  implicit none
  type(monolis_structure) :: A !> Ax = b の係数行列
  call monolis_set_tolerance(A, 1.0d-8) !> 反復法の収束判定閾値を 1.0e-8 回に設定
```
