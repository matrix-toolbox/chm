#!/usr/bin/env bash
# =============================================================================
# update_browse.sh -- rebuild matrices.js, the data behind browse.html
#
# 2026-09-13  Claude Opus 5
#
# matrices.js is generated, never edited by hand.  Run this after adding a
# matrix anywhere the Catalog records one:
#
#     index.html                 a new entry in the main list
#     CHM_dL/index.html          Appendix A: a probe, a family, a table row
#     CHM_SH/index.html          Appendix C: a symmetric or Hermitian row
#     CHM_*/*.m                  a new file named  <prefix>_N_d_L.m
#
#     ./update_browse.sh         rebuild and report what changed
#     ./update_browse.sh -n      only check whether it is out of date
#
# Exit status 1 in -n mode means matrices.js no longer matches the sources.
# =============================================================================
set -u
export LC_ALL=C
cd "$(dirname "$0")" || exit 1

CHECK=0
[ "${1:-}" = "-n" ] && CHECK=1

OUT=matrices.js
TMP=$(mktemp -t matrices.XXXXXX.js) || exit 1
trap 'rm -f "$TMP"' EXIT

STATS=$(python3 build_matrices.py -o "$TMP" | grep -v "^  wrote ") || exit 1

# the first line carries the build date, which changes every day -- compare
# only the data itself
same=0
if [ -f "$OUT" ] && [ "$(tail -n +2 "$OUT")" = "$(tail -n +2 "$TMP")" ]; then
    same=1
fi

if [ "$CHECK" = 0 ]; then
    echo "$STATS"
fi

if [ "$CHECK" = 1 ]; then
    if [ "$same" = 1 ]; then
        echo "matrices.js is up to date."
        exit 0
    fi
    echo "matrices.js is OUT OF DATE -- run ./update_browse.sh"
    exit 1
fi

if [ "$same" = 1 ]; then
    echo "no change: matrices.js already matches the sources."
    exit 0
fi

# what appeared and what went away, by name
names() { sed 's/^var matrices = //; s/;$//' "$1" |
          tr '{' '\n' | sed -n 's/.*"nm":"\([^"]*\)".*/\1/p' | sort -u; }
if [ -f "$OUT" ]; then
    added=$(comm -13 <(names "$OUT") <(names "$TMP") | tr '\n' ' ')
    gone=$(comm  -23 <(names "$OUT") <(names "$TMP") | tr '\n' ' ')
    [ -n "$added" ] && echo "  new names    : $added"
    [ -n "$gone"  ] && echo "  names gone   : $gone"
fi

cp "$TMP" "$OUT" || exit 1
echo "  updated $OUT -- reload browse.html"
