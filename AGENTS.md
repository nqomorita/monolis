# AGENTS.md

monolis (Metagraph-Oriented Network for Linear Iterative Solvers) は、MPI 並列の疎行列反復解法ライブラリです。
Fortran で実装され、C ラッパが提供されます。

## ビルドとテスト

サブモジュールと外部ライブラリ (METIS / ParMETIS / MUMPS / scalapack / monolis_utils / gedatsu / ggtools) の取得・ビルドが先に必要です。

```sh
./install_lib.sh        # 通常 (gfortran/OpenMPI)
./install_lib_intel.sh  # Intel oneAPI
./install_lib_A64FX.sh  # 富岳 (mpifrtpx)
make                    # ライブラリ + 全テストバイナリ
make FLAGS=DEBUG        # デバッグビルド (bounds check, fpe trap)
make FLAGS=INTEL        # Intel コンパイラビルド
make clean
```

ビルド成果物:
- `lib/libmonolis_solver.a` (Fortran 本体)、`lib/libmonolis.a` 相当のリンクは `USE_LIB1` を参照
- C ヘッダは `wrapper/*.h` から `include/` にコピーされる ([Makefile](Makefile) の `cp_header`)

テスト実行:
```sh
cd src_test && ./run.sh        # ローカル (Fortran + C)
cd src_test && ./run.CI.sh     # CI 用 (set -e、--oversubscribe)
```

CI のテストは [src_test/run.CI.sh](src_test/run.CI.sh) と [wrapper_test/run.CI.sh](wrapper_test/run.CI.sh)。
テスト関数の網羅性は [src_test/cover_check.sh](src_test/cover_check.sh) で `src/**/*.f90` の全 subroutine と実行ログを diff して確認します。

## ディレクトリ構成

| パス | 役割 |
|---|---|
| [src/](src/) | Fortran 本体。サブカテゴリは下表参照 |
| [wrapper/](wrapper/) | C API。`monolis_*_c.{c,h}` (C 側) と `*_wrapper.f90` (Fortran 側 `bind(c)` ブリッジ) のペア |
| [src_test/](src_test/) | Fortran ユニットテスト。`src/` と同じディレクトリ構造をミラー |
| [wrapper_test/](wrapper_test/) | C ユニットテスト |
| [submodule/](submodule/) | 外部依存。直接編集しない |
| [sample/](sample/) | 利用サンプル (mesher + solver) |
| [include/](include/) `lib/` `bin/` | ビルド出力。手で編集しない |

`src/` のサブカテゴリ:

| パス | 役割 |
|---|---|
| [src/define/](src/define/) | 中核データ構造の定義。パラメータ構造体 `monolis_prm` ([def_solver_prm.f90](src/define/def_solver_prm.f90))、行列構造体 `monolis_mat` ([def_mat.f90](src/define/def_mat.f90))、ソルバ全体構造体 `monolis_structure` ([def_struc.f90](src/define/def_struc.f90))、`Iarray`/`Rarray` のキー定数群 |
| [src/matrix/](src/matrix/) | 疎行列 (CSR) の組み立て・コピー・対称化・非ゼロパターン構築 (`spmat_nzpattern*`, `spmat_handler*`, `spmat_convert_sym`) |
| [src/linalg/](src/linalg/) | BLAS 風基本演算: 内積 (`inner_product`)、ベクトル操作 (`vec_util`)、行列ベクトル積 (`matvec`)、行列積 (`matmat`)、収束判定 (`mat_converge`)。MPI 通信込み |
| [src/prec/](src/prec/) | 前処理 (preconditioner)。`diag/` (対角)・`sor/` (SOR) のブロックサイズ別実装 (`*_33` 3x3、`*_nn` 一般、`*_V` 可変ブロック) と、ディスパッチ用 [precond.f90](src/prec/precond.f90) |
| [src/iterative/](src/iterative/) | 反復解法本体: CG / BiCGSTAB / GropCG / PipeCG / PipeCR / PipeBiCGSTAB / BiCGSAFE / IDRS / COCG / SOR。各ファイル 1 メソッドで `monolis_solver_<name>` を公開 |
| [src/solver/](src/solver/) | ソルバのトップエントリ [solver.f90](src/solver/solver.f90)。`monoPRM` のメソッド種別を見て `iterative/` のいずれかにディスパッチ |
| [src/eigen/](src/eigen/) | 固有値解法 (Lanczos 法とそのユーティリティ、`eigen_solver` のディスパッチ) |
| [src/optimize/](src/optimize/) | 最適化ルーチン。現状は非負最小二乗 [nnls.f90](src/optimize/nnls.f90) |
| [src/wrapper/](src/wrapper/) | 外部数値ライブラリ呼び出し (LAPACK / ScaLAPACK)。C ラッパの `wrapper/` ではなく、Fortran から外部 BLAS 系を叩くアダプタ |

新規ソースは `Makefile` の `SRC_*` リスト (および C 側 `SRC_*_C` / `SRC_*_CF`) と、対応するテストでは `TST_SRC_*` 系リストへの追加が必要です。Makefile はワイルドカードを使わず明示的に列挙しています。

## コーディング規約

- **言語**: Fortran (free form、`-std=legacy` 許容)。インデント 2 スペース、`implicit none` 必須。
- **モジュール命名**: ファイル `src/iterative/CG.f90` → `module mod_monolis_solver_CG`、公開ルーチン `monolis_solver_CG`。プレフィックス `monolis_` / `mod_monolis_` を一貫使用。
- **doc コメント**: Doxygen 形式。`!>` で 1 行説明、引数は `!> [in]` `!> [in,out]` で記述。`@ingroup` でグループ分け。日本語コメントを使用。
- **型**: `kint` (int), `kdouble` (real(8)), C 側は `kint_c` / `c_double`。`mod_monolis_utils` 由来。
- **構造体引数の慣習**: ソルバルーチンは `(monoPRM, monoCOM, monoMAT, monoPREC)` の順で受け取る ([src/iterative/CG.f90](src/iterative/CG.f90) 参照)。
- **配列確保**: `monolis_alloc_R_1d` 等のユーティリティを使用 (ベア `allocate` ではなく)。
- **ベクトル/内積**: `monolis_inner_product_main_R`、`monolis_vec_AXPBY_R`、`monolis_matvec_product_main_R` 等の既存 BLAS 風ラッパを使用。
- **デバッグログ**: 各ソルバ冒頭で `call monolis_std_debug_log_header("...")`。

## C ラッパパターン

新しい Fortran ルーチンを C から呼べるようにする際は:
1. `wrapper/<area>/<name>_wrapper.f90` に `bind(c, name="monolis_<name>_c_main")` ルーチンを追加。
2. `wrapper/<area>/monolis_<name>_c.c` でユーザ向け C 関数を定義し、`monolis_<name>_c_main` を呼ぶ。
3. ヘッダ `monolis_<name>_c.h` を作り `wrapper/monolis_solver.h` (相当) に取り込む。
4. [Makefile](Makefile) の `SRC_*_CF` (Fortran) と `SRC_*_C` (C) リストに追加。
5. テストは [wrapper_test/](wrapper_test/) に同名で追加し `SRC_*_C_TEST` に登録。

C 側は **0-based**、Fortran 側に渡す際に通信テーブル (`recv_item`/`send_item`) は `+1` で 1-based に変換します ([wrapper/solver/solver_wrapper.f90](wrapper/solver/solver_wrapper.f90))。

## テスト追加のルール

- 新規 subroutine を `src/` に追加したら、必ず `src_test/<同パス>/<name>_test.f90` を追加し、内部で `monolis_std_global_log_string("<関数名>")` を呼ぶ。`cover_check.sh` がこのログ出力で網羅率を判定します。
- テストモジュール命名: `mod_monolis_<name>_test`、エントリ `monolis_<name>_test`、これを [src_test/monolis_test.f90](src_test/monolis_test.f90) から呼び出します。
- `n_dof = 1..3` と前処理 `monolis_prec_NONE/DIAG/SOR` の組合せ走査が反復解法テストの定型 ([src_test/iterative/CG_test.f90](src_test/iterative/CG_test.f90))。

## 注意点

- **submodule は編集しない**。`monolis_utils` / `gedatsu` の API 変更が必要な場合はそのリポジトリで作業。
- `install_lib.sh` の再実行は scalapack / MUMPS の cmake 再構築が走り時間がかかる。サブモジュールのみ再ビルドしたければ該当ディレクトリで `make` を直接実行。
- ファイル `monolis_solver.f90` と `monolis.f90` は集約モジュール。新規 use 文の追加先はここ。
- Doxygen マニュアルは `manual/{c,fortran}/doxyfile` から生成。コメントの doxygen 形式を崩さないこと。
