# monolis

- Morita's non-overlapping / overlapping DDM based linear equation solver

## インストール

### 0. インストールの準備

Gitlab の [SSH keys](https://gitlab.com/profile/keys) から、公開鍵を登録する。
例えば、[Qiita:【GitLab】SSH認証キー（公開鍵）を登録する](https://qiita.com/CUTBOSS/items/462a2ed28d264aeff7d5) などに詳しい。

### 1. クローン

Gitlab から monolis リポジトリをローカルにクローンする。

```
$ git clone git@gitlab.com:morita/monolis.git
```

monolis ディレクトリに移動する。

```
$ cd monolis
```

### 2. インストール

- (A) 逐次計算で利用する場合

make する。

```
$ make
```

- (B) グラフライブラリ metis を利用する場合


Makefile の中の以下の部分を適切に設定する。

```
METIS_DIR = $(path to metis)
METIS_INC = -I $(METIS_DIR)/include
METIS_LIB = -L$(METIS_DIR)/lib -lmetis
```

make する。


```
$ make FLAGS=METIS
```

- (C) 並列計算で利用する場合

Makefile の中の以下の部分を適切に設定する。

```
FC     = mpif90
```

make する。

```
$ make FLAGS=MPI
```

- (D) 全てを利用する場合

make する。

```
$ make FLAGS=MPI,METIS
```

### 3. インストールの確認

インストールが成功していれば、`monolis/lib` に `libmonolis.a` が生成されている。
以下のコマンドで確認できる。

```
$ ls lib
```

