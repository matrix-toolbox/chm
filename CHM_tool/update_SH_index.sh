#!/usr/bin/env bash
# =============================================================================
# first: ./show_histogram.sh 
# update_SH_index.sh -- merge dL_statistic.md into Appendix C (CHM_SH/index.html)
#
# 2026-09-11  Claude Opus 5
#
# Reads the (N, symmetry, d, #Lambda) table produced by ./show_histogram.sh
# and adds, to the matching table of CHM_SH/index.html, every (d, #Lambda)
# pair that is not listed there yet.  Added rows are greyed out (#aaa), have
# empty "matrix" and "BH class" cells, and carry the comment
#
#     observed numerically
#
# Nothing existing is ever reordered or removed -- new rows are merged into
# place so each table stays sorted by d, then #Lambda.
#
#   symmetry = s  ->  the "Symmetric CHM" section
#   symmetry = h  ->  the "Hermitian CHM" section
#
# The two catch-all tables ("N > 11" for symmetric, "N > 12" for Hermitian)
# are SPLIT into one table per dimension.  The N of each existing row is
# taken from its matrix label (F<sub>13</sub> -> 13, RH<sub>20,171,2</sub>
# -> 20); rows that carry no such label (the "..." placeholders) are kept in
# a final "N > <largest>" table.  A dimension seen only in dL_statistic.md
# gets a table of its own in the right place, so every row lives under an
# explicit N and no comment has to spell the dimension out.
#
# Running this again on an already-split page is a no-op structurally.
#
# A few HTML errors in the current page are repaired on the way (see --help).
#
#   ./update_SH_index.sh                       # update ../CHM_SH/index.html
#   ./update_SH_index.sh -o /tmp/preview.html  # leave the page alone
#
# The page is rewritten in place.  The script is idempotent -- a second run on
# an already-updated page changes nothing -- and the page is under version
# control, so `git diff CHM_SH/index.html` shows exactly what a run did.
# =============================================================================

set -u -o pipefail
export LC_ALL=C

STATS="dL_statistic.md"
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SRC="$HERE/../CHM_SH/index.html"
OUT=""                  # empty = rewrite SRC in place
LIVE_URL="https://matrix-toolbox.github.io/chm/CHM_SH/index.html"

usage() {
    cat <<'EOF'
usage: update_SH_index.sh [-s STATS] [-i INDEX] [-o OUT]

  -s STATS   markdown written by ./show_histogram.sh   (default dL_statistic.md)
  -i INDEX   the Appendix C page to read       (default ../CHM_SH/index.html;
             if that file does not exist the published page is downloaded)
  -o OUT     file to write                     (default: rewrite INDEX itself)

HTML repairs always applied:
  * "<p></p>Once a matrix ... </p>"  ->  "<p>Once a matrix ... </p>"   (stray </p>)
  * "<tr><td><td>"                   ->  "<tr><td></td><td>"           (unclosed <td>)
  * <h1>Catalog of CHM</h1>          ->  <h1 id="top">...              (the
    "go top" arrows all point at #top, which did not exist)
  * <a href="">SH<sub>7,0,97</sub></a> and the like get the href they
    obviously meant (derived from the label)
  * a closing </html> is appended when missing
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        -s) STATS="$2"; shift 2 ;;
        -i) SRC="$2";   shift 2 ;;
        -o) OUT="$2";   shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "update_SH_index.sh: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

[ -f "$STATS" ] || { echo "update_SH_index.sh: no such file: $STATS" >&2; exit 2; }

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# no local copy of the page?  take the published one, so the script works
# straight out of the toolbox directory.
if [ ! -f "$SRC" ]; then
    echo "update_SH_index.sh: $SRC not found, fetching the published page" >&2
    if ! curl -fsS -m 30 "$LIVE_URL" -o "$TMP/live.html"; then
        echo "update_SH_index.sh: could not fetch $LIVE_URL (pass -i FILE instead)" >&2
        exit 2
    fi
    SRC="$TMP/live.html"
    [ -n "$OUT" ] || OUT="index_NEW.html"   # nothing local to rewrite
fi
[ -n "$OUT" ] || OUT="$SRC"

# ---- 1. repair the known HTML faults, before anything else ----------------
awk '
function fix_empty_href(line,   out, pre, m, label, fname) {
    out = ""
    while (match(line, /<a href="">[A-Za-z]+<sub>[^<]*<\/sub><\/a>/)) {
        pre  = substr(line, 1, RSTART - 1)
        m    = substr(line, RSTART, RLENGTH)
        line = substr(line, RSTART + RLENGTH)
        label = m
        sub(/^<a href="">/, "", label); sub(/<\/a>$/, "", label)
        fname = label
        gsub(/<sub>/, "_", fname); gsub(/<\/sub>/, "", fname); gsub(/,/, "_", fname)
        out = out pre "<a href=\"https://github.com/matrix-toolbox/CHM_CATALOG/blob/master/CHM_SH/" fname ".m\">" label "</a>"
    }
    return out line
}
{
    line = $0
    sub(/<p><\/p>Once a matrix/, "<p>Once a matrix", line)      # stray </p>
    gsub(/<tr><td><td>/, "<tr><td></td><td>", line)             # unclosed <td>
    sub(/<h1>Catalog of CHM<\/h1>/, "<h1 id=\"top\">Catalog of CHM</h1>", line)
    if (line ~ /<a href="">/) line = fix_empty_href(line)
    print line
}
' "$SRC" > "$TMP/fixed.html"

# add the grey style for the generated rows, next to the existing .comment rule
if ! grep -q '^\.observed{' "$TMP/fixed.html"; then
    sed -i 's#^\.comment{color:\#aaa;}#.comment{color:\#aaa;}\n.observed{color:\#aaa;}#' "$TMP/fixed.html"
fi

# ---- 2. what does the page already contain? ------------------------------
#   sections.txt : "sec N"        every <h3>N = k</h3> table
#   have.txt     : "sec N d L"    every listed pair; for catch-all rows the N
#                                 is taken from the matrix label
awk '
function label_N(cell,   s) {                 # F<sub>13</sub> -> 13
    s = cell
    if (match(s, /<sub>[0-9]+/)) {
        s = substr(s, RSTART + 5, RLENGTH - 5)
        return s
    }
    return ""
}
/id="Symmetric_CHM"/ {sec="s"}
/id="Hermitian_CHM"/ {sec="h"}
/<h3[^>]*>N = [0-9]+<\/h3>/ { t=$0; sub(/.*N = /,"",t); sub(/<\/h3>.*/,"",t); N=t; catchall=0; print sec, N > SECT }
/<h3[^>]*>N &gt;/ { N=""; catchall=1 }
/<tbody>/ {inb=1; next}
/<\/tbody>/ {inb=0; next}
inb && /<tr/ {
    n=split($0, p, /<td[^>]*>/)
    c1=p[2]; sub(/<\/td>.*/,"",c1)
    d=p[3]; sub(/<\/td>.*/,"",d); gsub(/<[^>]*>/,"",d); gsub(/[ \t]/,"",d)
    L=p[4]; sub(/<\/td>.*/,"",L); gsub(/<[^>]*>/,"",L); gsub(/[ \t]/,"",L)
    rowN = catchall ? label_N(c1) : N
    if (rowN != "" && d ~ /^[0-9]+$/ && L ~ /^[0-9]+$/) print sec, rowN, d, L > HAVE
}
' SECT="$TMP/sections.txt" HAVE="$TMP/have.txt" "$TMP/fixed.html"
touch "$TMP/sections.txt" "$TMP/have.txt"

# ---- 3. what does dL_statistic.md claim?   want.txt : "sec N d L" ---------
awk '
/^## N = / {
    t=$0
    sub(/^## N = /,"",t); sym=t; sub(/,.*/,"",t); N=t
    sub(/.*symmetry = /,"",sym); sub(/[^a-z].*/,"",sym)
    next
}
/^\| *[0-9]+ *\| *[0-9]+ *\|/ {
    split($0, c, "|")
    d=c[2]; L=c[3]; gsub(/[ \t]/,"",d); gsub(/[ \t]/,"",L)
    if (N != "" && sym != "") print sym, N, d, L
}
' "$STATS" > "$TMP/want.txt"

# ---- 4. which pairs are missing?   insert.txt : "sec N d L" --------------
awk '
FILENAME==HAVE { listed[$1 SUBSEP $2 SUBSEP $3 SUBSEP $4]=1; next }
{ if (!(($1 SUBSEP $2 SUBSEP $3 SUBSEP $4) in listed)) print $1, $2, $3, $4 }
' HAVE="$TMP/have.txt" "$TMP/have.txt" "$TMP/want.txt" \
  | sort -k1,1 -k2,2n -k3,3n -k4,4n > "$TMP/insert.txt"

N_NEW=$(wc -l < "$TMP/insert.txt")

# ---- 5. write the new page ------------------------------------------------
awk -v INS="$TMP/insert.txt" -v SECT="$TMP/sections.txt" '
function key_gt(ad, aL, bd, bL) {             # is (ad,aL) > (bd,bL) ?
    if (bd == "") return 0                    # non-numeric row: keep new first
    if (ad + 0 != bd + 0) return (ad + 0) > (bd + 0)
    return (aL + 0) > (bL + 0)
}
function new_row(d, L) {
    return "<tr class=\"observed\"><td></td><td>" d "</td><td>" L "</td><td></td><td>observed numerically</td></tr>"
}
function label_N(cell,   s) {
    s = cell
    if (match(s, /<sub>[0-9]+/)) return substr(s, RSTART + 5, RLENGTH - 5)
    return ""
}
function row_d(line,   p, n, d) {
    n = split(line, p, /<td[^>]*>/); d = p[3]
    sub(/<\/td>.*/, "", d); gsub(/<[^>]*>/, "", d); gsub(/[ \t]/, "", d)
    return (d ~ /^[0-9]+$/) ? d : ""
}
function row_L(line,   p, n, L) {
    n = split(line, p, /<td[^>]*>/); L = p[4]
    sub(/<\/td>.*/, "", L); gsub(/<[^>]*>/, "", L); gsub(/[ \t]/, "", L)
    return (L ~ /^[0-9]+$/) ? L : ""
}
function row_c1(line,   p, n, c) {
    n = split(line, p, /<td[^>]*>/); c = p[2]
    sub(/<\/td>.*/, "", c)
    return c
}
# merge the buffered rows of one table with the pending inserts for (sec,N)
function merge(s, nn, cnt_rows, txt, dd, ll,    i, j, k) {
    k = s SUBSEP nn
    j = 1
    for (i = 1; i <= cnt_rows; i++) {
        while (j <= icnt[k] && !key_gt(ins_d[k,j], ins_L[k,j], dd[i], ll[i])) {
            print new_row(ins_d[k,j], ins_L[k,j]); j++
        }
        print txt[i]
    }
    while (j <= icnt[k]) { print new_row(ins_d[k,j], ins_L[k,j]); j++ }
}
function emit_split(   i, j, n, nlist, cnt, tmp, maxN, k, ii) {
    # dimensions to emit: those found in the catch-all rows, plus any that
    # dL_statistic.md mentions and that has no <h3>N = k</h3> of its own
    cnt = 0
    for (n in split_has) { if (split_has[n] == sec) { cnt++; nlist[cnt] = n } }
    for (k in icnt) {
        split(k, kk, SUBSEP)
        if (kk[1] != sec) continue
        if ((sec SUBSEP kk[2]) in dedicated) continue
        already = 0
        for (ii = 1; ii <= cnt; ii++) if (nlist[ii] == kk[2]) already = 1
        if (!already) { cnt++; nlist[cnt] = kk[2] }
    }
    for (i = 2; i <= cnt; i++) {                       # insertion sort, ascending
        tmp = nlist[i]
        for (j = i - 1; j >= 1 && nlist[j] + 0 > tmp + 0; j--) nlist[j+1] = nlist[j]
        nlist[j+1] = tmp
    }
    # the tail heading must never go backwards: on a page that is already
    # split there is nothing left to split out, and the threshold has to stay
    # whatever it was ("N > 20"), not collapse to "N > 0".
    maxN = orig_thresh + 0
    if (cnt > 0 && nlist[cnt] + 0 > maxN) maxN = nlist[cnt] + 0
    for (i = 1; i <= cnt; i++) {
        n = nlist[i]
        print "<h3>N = " n "</h3>"
        print "<table>"
        for (j = 1; j <= nthead; j++) print thead[j]
        print "<tbody>"
        ncnt = srow_cnt[sec SUBSEP n]
        for (j = 1; j <= ncnt; j++) {
            mtxt[j] = srow_txt[sec SUBSEP n, j]
            md[j]   = srow_d[sec SUBSEP n, j]
            ml[j]   = srow_l[sec SUBSEP n, j]
        }
        merge(sec, n, ncnt, mtxt, md, ml)
        print "</tbody>"
        print "</table>"
        print ""
    }
    if (ntail > 0) {
        print "<h3>N &gt; " maxN "</h3>"
        print "<table>"
        for (j = 1; j <= nthead; j++) print thead[j]
        print "<tbody>"
        for (j = 1; j <= ntail; j++) print tail[j]
        print "</tbody>"
        print "</table>"
        # no trailing blank line here: whatever followed </table> in the
        # source still prints, so the spacing stays as it was (and does not
        # grow by one blank line on every re-run)
    }
}
BEGIN {
    while ((getline line < INS) > 0) {
        split(line, a, " ")
        k = a[1] SUBSEP a[2]
        icnt[k]++
        ins_d[k, icnt[k]] = a[3]; ins_L[k, icnt[k]] = a[4]
    }
    close(INS)
    while ((getline line < SECT) > 0) { split(line, a, " "); dedicated[a[1] SUBSEP a[2]] = 1 }
    close(SECT)
}
/id="Symmetric_CHM"/ { sec = "s" }
/id="Hermitian_CHM"/ { sec = "h" }

# ---- the catch-all table: capture it whole, then re-emit it split by N ----
!catch_mode && /<h3[^>]*>N &gt;/ {
    catch_mode = 1; nthead = 0; ntail = 0; in_thead = 0; in_tb = 0
    orig_thresh = $0                                  # remember "N &gt; 12"
    sub(/.*N &gt;[ ]*/, "", orig_thresh); sub(/<\/h3>.*/, "", orig_thresh)
    next
}
catch_mode {
    if ($0 ~ /<thead>/)  { in_thead = 1 }
    if (in_thead)        { thead[++nthead] = $0 }
    if ($0 ~ /<\/thead>/){ in_thead = 0; next }
    if (in_thead)        { next }
    if ($0 ~ /<tbody>/)  { in_tb = 1; next }
    if ($0 ~ /<\/tbody>/){ in_tb = 0; next }
    if (in_tb && $0 ~ /<tr/) {
        n = label_N(row_c1($0))
        if (n != "") {
            k = sec SUBSEP n
            srow_cnt[k]++
            srow_txt[k, srow_cnt[k]] = $0
            srow_d[k, srow_cnt[k]]   = row_d($0)
            srow_l[k, srow_cnt[k]]   = row_L($0)
            split_has[n] = sec
        } else {
            dup = 0
            for (i = 1; i <= ntail; i++) if (tail[i] == $0) dup = 1
            if (!dup) tail[++ntail] = $0
        }
        next
    }
    if ($0 ~ /<\/table>/) { emit_split(); catch_mode = 0 }
    next
}

# ---- ordinary <h3>N = k</h3> tables: merge in place -----------------------
/<h3[^>]*>N = [0-9]+<\/h3>/ { t = $0; sub(/.*N = /, "", t); sub(/<\/h3>.*/, "", t); N = t }
/<tbody>/  { print; inb = 1; nbuf = 0; next }
/<\/tbody>/ {
    for (i = 1; i <= nbuf; i++) { btxt[i] = buf[i]; bd[i] = bufd[i]; bl[i] = bufl[i] }
    merge(sec, N, nbuf, btxt, bd, bl)
    print; inb = 0; next
}
inb {
    nbuf++; buf[nbuf] = $0
    bufd[nbuf] = ($0 ~ /<tr/) ? row_d($0) : ""
    bufl[nbuf] = ($0 ~ /<tr/) ? row_L($0) : ""
    next
}
{ print }
' "$TMP/fixed.html" > "$TMP/out.html"

grep -qi '</html>' "$TMP/out.html" || printf '</html>\n' >> "$TMP/out.html"

# built aside and moved into place, so an interrupted run cannot truncate the page
mv -- "$TMP/out.html" "$OUT"

echo "wrote $OUT  (+${N_NEW} rows from $STATS)"
if [ "$N_NEW" -gt 0 ]; then
    echo "added:"
    awk '{printf "  symmetry=%s  N=%-4s d=%-4s #L=%s\n", $1, $2, $3, $4}' "$TMP/insert.txt"
fi
