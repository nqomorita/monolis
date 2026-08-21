#!/bin/bash
#> Fortran の module / use 関係を解析して make の依存関係 (.depend) を生成する
#> 並列ビルド (make -j) で .mod の生成順序を保証するために使用する

set -u

#> ソース中の日本語コメントで awk がマルチバイト解釈エラーを起こすため C ロケールで走らせる
export LC_ALL=C

DIRS="src wrapper src_test driver"

FILES=""
for d in $DIRS; do
  if [ -d "$d" ]; then
    FILES="$FILES $(find "$d" -name '*.f90')"
  fi
done

if [ -z "$FILES" ]; then
  exit 0
fi

awk '
#> ソースパスの先頭ディレクトリを obj に置き換えてオブジェクトパスを得る
function objof(path,   n, a, i, r) {
  n = split(path, a, "/")
  r = "obj"
  for (i = 2; i <= n; i++) r = r "/" a[i]
  sub(/\.f90$/, ".o", r)
  return r
}
FNR == 1 { obj = objof(FILENAME); objs[FILENAME] = obj }
{
  line = tolower($0)
  if (match(line, /^[ \t]*module[ \t]+[a-z][a-z0-9_]*/)) {
    name = substr(line, RSTART, RLENGTH)
    sub(/^[ \t]*module[ \t]+/, "", name)
    if (name != "procedure") provider[name] = obj
  }
  else if (match(line, /^[ \t]*use[ \t]+(::[ \t]*)?[a-z][a-z0-9_]*/)) {
    name = substr(line, RSTART, RLENGTH)
    sub(/^[ \t]*use[ \t]+(::[ \t]*)?/, "", name)
    uses[obj] = uses[obj] " " name
  }
}
END {
  for (o in uses) {
    n = split(uses[o], u, " ")
    delete seen
    dep = ""
    for (i = 1; i <= n; i++) {
      m = u[i]
      if (m == "" || !(m in provider)) continue
      p = provider[m]
      if (p == o || (p in seen)) continue
      seen[p] = 1
      dep = dep " " p
    }
    if (dep != "") print o ":" dep
  }
}
' $FILES | sort
