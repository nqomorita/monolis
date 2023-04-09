#!/bin/bash

# テストで実行された関数名のリストを作成
grep MONOLIS test_list.dat | sort | uniq | sed "s/\*\* MONOLIS: //g" | grep -v \\[ | grep -v "** MONOLIS " | grep -v \\[ | grep -v "DOF  " | grep -v "PRECOND  " > 1_uniq_test_list.dat

# monolis で実装された関数名のリストを作成
cd ..

find ./src -name "*.f90" | xargs grep -ih subroutine | grep -v "  end" | grep -v \! | sed "s/  subroutine //g" | sed "s/  recursive subroutine //g" | sed "s/(.*//g" > ./src_test/tmp_list.dat

find ./src -name "*.f90" | xargs grep -ih function | grep -v end | grep -v \! | sed "s/  function //g" | sed "s/(.*//g" >> ./src_test/tmp_list.dat

cd ./src_test

sort tmp_list.dat > 0_impl_list.dat

# 差分チェック
diff -i 0_impl_list.dat 1_uniq_test_list.dat > 2_diff
