# monolis

- Morita's non-overlapping / overlapping DDM based linear equation solver

<img src="./monolis.logo.svg" height=100px>

## インストール

### 0. インストールの準備

Gitlab の [SSH keys](https://gitlab.com/profile/keys) から、公開鍵を登録する。
例えば、[Qiita:【GitLab】SSH認証キー（公開鍵）を登録する](https://qiita.com/CUTBOSS/items/462a2ed28d264aeff7d5) などに詳しい。

### 1. インストール環境の準備

インストールには以下の環境が必要である。

- make
- cmake
- git
- gcc (gfortran)
- MPI

Ubuntu 環境では、以下のコマンドでインストールできる。

```bash
sudo apt update
sudo apt upgrade
sudo apt install -y build-essential
sudo apt install -y cmake
sudo apt install -y gfortran
sudo apt install -y git
sudo apt install -y wget
sudo apt install -y openmpi-doc openmpi-bin libopenmpi-dev
```
CentOS 環境では、以下のコマンドでインストールできる。

```bash
sudo yum update
sudo yum groupinstall -y "Development Tools"
sudo yum install -y make
sudo yum install -y cmake
sudo yum install -y git
sudo yum install -y wget
sudo yum install -y openmpi openmpi-devel
```

### 2. クローン

Gitlab から monolis リポジトリをローカルにクローンする。

```bash
$ git clone git@gitlab.com:morita/monolis.git
```

monolis ディレクトリに移動する。

```bash
$ cd monolis
```

### 3. ライブラリのインストール

monolis にはグラフ分割ライブラリ metis を用いる。metis のインストールのために、monolis ディレクトリにおいて以下のコマンドを実行する。

```bash
$ ./install_lib.sh
```

### 4. monolis のインストール

<!--
- (A) 逐次計算で利用する場合

make する。

```bash
$ make
```
- (B) 並列計算で利用する場合

Makefile の中の以下の部分を適切に設定する。

```bash
FC     = mpif90
```
-->

monolis ディレクトリにおいて以下のコマンドを実行する。

```bash
$ make
```

### 3. インストールの確認

インストールが成功していれば、`monolis/bin` に `monolis_partitioner`、`monolis/lib` に `libmonolis.a` が生成されている。
以下のコマンドで確認できる。

```bash
$ ls bin
$ ls lib
```
