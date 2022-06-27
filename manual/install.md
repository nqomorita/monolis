
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

monolis ディレクトリにおいて以下のコマンドを実行する。

```bash
$ make
```

`make` と実行した場合は、MPIによる並列計算機能と、グラフ分割ライブラリMETISが有効になる。
`make` 実行時のオプションは以下の通り。

- MPI: MPI による並列計算の有効
- OPENMP: OpenMP による並列計算の有効
- METIS: グラフ分割ライブラリ METIS を利用
- MUMPS: 直接法ライブラリ MUMPS を利用
- DEBUG: DEBUG オプションを有効

実行例：MPI,OpenMPによる並列計算と、METIS、MUMPS ライブラリを有効にする場合

```bash
$ make FLAG=MPI,OPENMP,METIS,MUMPS
```


### 3. インストールの確認

インストールが成功していれば、`monolis/bin` に `monolis_partitioner`、`monolis/lib` に `libmonolis.a` が生成されている。
以下のコマンドで確認できる。

```bash
$ ls bin
$ ls lib
```


### 4. git による monolis のアップデート

monolis にアップデートが生じた場合は、monolis トップディレクトリにおいて、
以下の git コマンドを利用することにより最新のソースにアップデートできる。

```
$ make clean
$ git pull
$ make
```
