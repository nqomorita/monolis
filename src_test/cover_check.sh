#!/bin/bash

# テストで実行された関数名のリストを作成
grep MONOLIS test_list.dat | sort | uniq | sed "s/\*\* MONOLIS: //g" | grep -v \\[ | grep -v "** MONOLIS " | grep -v \\[ | grep -v "DOF  " | grep -v "PRECOND  " > 1_uniq_test_list.dat

# monolis で実装された関数名のリストを作成
# モジュールレベル (2 スペースインデント) の subroutine のみ対象とする
# (interface 宣言や contains 内部手続きは 4 スペース以上のため除外される)
cd ..

find ./src -name "*.f90" | xargs grep -hE "^  (recursive )?subroutine [a-zA-Z0-9_]+" | sed -E "s/^  (recursive )?subroutine +([a-zA-Z0-9_]+).*/\2/" > ./src_test/tmp_list.dat

find ./src -name "*.f90" | xargs grep -hiE "^  [^!]*function +[a-zA-Z0-9_]+" | grep -iv "end function" | sed -E "s/.*function +([a-zA-Z0-9_]+).*/\1/" >> ./src_test/tmp_list.dat

cd ./src_test

sort tmp_list.dat > 0_impl_list.dat

# 差分チェック
diff -i 0_impl_list.dat 1_uniq_test_list.dat > 2_diff
