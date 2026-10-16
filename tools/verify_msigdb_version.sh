#!/usr/bin/env bash
# =====================================================================
# verify_msigdb_version.sh — 用字节级比对**定版**本地 MSigDB gmt
#
# 背景：gmt 文件内部不含版本标记，早期只能写「不写版本 + 报规模与校验和」。
#       但 MSigDB 的历史归档在 Broad 公开目录下**无需登录**即可取到：
#         https://data.broadinstitute.org/gsea-msigdb/msigdb/release/<版本>/<集合>.v<版本>.symbols.gmt
#       因此可以用 md5 把本地文件与各版本**逐字节**比对，从而确证版本。
#
# 用法： bash tools/verify_msigdb_version.sh
# 输出： data/_vercheck/msigdb_version_verdict.txt（同时打印到 stdout）
# 只读远端 + 只写 _vercheck 缓存，不改动 data/genesets/ 下的正式输入。
# =====================================================================
set -u
cd "$(dirname "$0")/.." || exit 1
BASE="https://data.broadinstitute.org/gsea-msigdb/msigdb/release"
CACHE="data/_vercheck"
OUT="$CACHE/msigdb_version_verdict.txt"
VERSIONS="7.5.1 2023.1.Hs 2023.2.Hs 2024.1.Hs 2025.1.Hs 2026.1.Hs"
SETS="c5.go.bp h.all"

mkdir -p "$CACHE"
{
  echo "MSigDB 版本判定报告  生成时间：$(date '+%Y-%m-%d %H:%M:%S %z')"
  echo "方法：本地 data/genesets/*.gmt 与官方公开归档的文件做 md5 逐字节比对"
  echo "归档：$BASE/<版本>/<集合>.v<版本>.symbols.gmt"
  echo

  echo "== 本地文件指纹 =="
  for s in $SETS; do
    f="data/genesets/$s.symbols.gmt"
    if [ -f "$f" ]; then
      printf "  %-28s %10d B  %s  sets=%s\n" "$s.symbols.gmt" "$(stat -c%s "$f")" \
        "$(md5sum "$f" | cut -d' ' -f1)" "$(wc -l < "$f")"
    else
      printf "  %-28s 缺失\n" "$s.symbols.gmt"
    fi
  done
  echo

  echo "== 官方归档各版本指纹 =="
  printf "  %-11s %-10s %-34s %10s %6s\n" "版本" "集合" "md5" "字节" "sets"
  for v in $VERSIONS; do
    for s in $SETS; do
      f="$CACHE/r_$s.v$v.symbols.gmt"
      if [ ! -f "$f" ]; then
        curl -s --noproxy '*' --max-time 300 -o "$f" "$BASE/$v/$s.v$v.symbols.gmt" || true
      fi
      if [ -s "$f" ] && [ "$(stat -c%s "$f")" -gt 10000 ]; then
        printf "  %-11s %-10s %-34s %10d %6d\n" "$v" "$s" \
          "$(md5sum "$f" | cut -d' ' -f1)" "$(stat -c%s "$f")" "$(wc -l < "$f")"
      else
        printf "  %-11s %-10s %s\n" "$v" "$s" "（归档不可取或为空）"
      fi
    done
  done
  echo

  echo "== 判定 =="
  for s in $SETS; do
    lf="data/genesets/$s.symbols.gmt"
    [ -f "$lf" ] || continue
    lm=$(md5sum "$lf" | cut -d' ' -f1)
    hit=""
    n=0
    for v in $VERSIONS; do
      f="$CACHE/r_$s.v$v.symbols.gmt"
      [ -s "$f" ] || continue
      n=$((n+1))
      if [ "$(md5sum "$f" | cut -d' ' -f1)" = "$lm" ]; then hit="$hit $v"; fi
    done
    if [ -n "$hit" ]; then
      echo "  $s.symbols.gmt  →  命中版本：$hit   （在 $n 个候选版本中比对）"
    else
      echo "  $s.symbols.gmt  →  未命中任何候选版本（md5 $lm）"
    fi
  done
} | tee "$OUT"

echo
echo "→ 报告已写入 $OUT"
