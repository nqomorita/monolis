
## 既存ソルバへの組み込み

monolis の関数を使う場合は、`mod_monolis` モジュールを `use` することで利用することができる。

monolis を使った連立一次方程式の求解はおもに以下の手順で実現される。

1. monolis 全体の初期化処理：[initialize](./impl/init.md)
1. monolis 構造体の初期化処理：[initialize](./impl/init.md)
1. 疎行列の非零パターンの決定：[matrix init](./impl/matrix_init.md)
1. 全体剛性行列への値の足し込み：[matrix add](./impl/matrix_add.md)
1. 境界条件の追加：[solver](./impl/solve.md)
1. 線形ソルバのパラメータ設定：[parameter](./impl/parameter.md)
1. 線形ソルバによる求解：[solver](./impl/solve.md)
1. monolis 構造体の終了処理：[initialize](./impl/init.md)
1. monolis 全体の修了処理：[initialize](./impl/init.md)


## 関数一覧

- [initialize](./impl/init.md)
- [matrix init](./impl/matrix_init.md)
- [matrix add](./impl/matrix_add.md)
- [solver](./impl/solve.md)
- [linear algebra](./impl/la.md)
- [parameter](./impl/parameter.md)
- [utilities](./impl/util.md)
