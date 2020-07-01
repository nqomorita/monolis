# monolis

- Morita's non-overlapping / overlapping DDM based linear equation solver

<img src="./monolis.logo.svg" height=100px>

## インストール

### 0. インストールの準備

Gitlab の [SSH keys](https://gitlab.com/profile/keys) から、公開鍵を登録する。
例えば、[Qiita:【GitLab】SSH認証キー（公開鍵）を登録する](https://qiita.com/CUTBOSS/items/462a2ed28d264aeff7d5) などに詳しい。

### 1. クローン

Gitlab から monolis リポジトリをローカルにクローンする。

```bash
$ git clone git@gitlab.com:morita/monolis.git
```

monolis ディレクトリに移動する。

```bash
$ cd monolis
```

### 2. ライブラリのインストール

monolis にはグラフ分割ライブラリ metis を用いる。metis のインストールのために、monolis ディレクトリにおいて以下のコマンドを実行する。

```bash
$ ./install_lib.sh
```

### 3. monolis のインストール

- (A) 逐次計算で利用する場合

make する。

```bash
$ make FLAGS=METIS
```

- (B) 並列計算で利用する場合

Makefile の中の以下の部分を適切に設定する。

```bash
FC     = mpif90
```

make する。

```bash
$ make FLAGS=MPI,METIS
```

### 3. インストールの確認

インストールが成功していれば、`monolis/bin` に `monolis_partitioner`、`monolis/lib` に `libmonolis.a` が生成されている。
以下のコマンドで確認できる。

```bash
$ ls bin
$ ls lib
```
