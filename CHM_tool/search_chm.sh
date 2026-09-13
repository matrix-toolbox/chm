#!/usr/bin/env bash
# =============================================================================
# search_chm.sh -- collect K complex Hadamard matrices whose Haagerup
#                  cardinality #L falls in a given range
#
# 2026-09-04  W. Bruzda  <w.bruzda@cft.edu.pl>
#
# Repeats
#     ./get_chm --size N --symmetry X
# and keeps only those matrices with  L1 <= #L <= L2,  until K of them have been
# collected (or an attempt / time limit is reached).  Rejected matrices are
# deleted unless --keep-rejects is given.
#
# Requires get_chm, lambda and (with --verify) summary and defect in the same
# directory, or on PATH via --bindir.
#
# 2026-09-13  every run collects its matrices into a single ../CHM_NUMERICAL,
#             instead of a chm_N<N>_<symmetry> directory per (size, symmetry).
#             The catalog, CHM_NUMERICAL.tsv, stays next to this script: the
#             matrices are disposable, the statistics are not.
#
#   ./search_chm.sh --count 10 --lmin 1 --lmax 100
#   ./search_chm.sh -k 5 -l 1500 -u 1600 --size 11 --symmetry s
#   ./search_chm.sh -k 3 -l 0 -u 200 --max-seconds 3600 --seed 20260904
# =============================================================================

set -u -o pipefail

# ---------------------------------------------------------------- defaults ---
COUNT=1                 # how many keepers to collect
LMIN=0                  # inclusive lower bound on #L
LMAX=-1                 # inclusive upper bound on #L; -1 = no upper bound
SIZE=11
SYMMETRY=s
METHOD=sinkhorn
ITERATIONS=1e5          # 1e+4 (the get_chm default) is thin above N = 12
LAMBDA_EPS=""           # empty = let lambda use its own default
OUTDIR=""               # empty = ../CHM_NUMERICAL
MAX_ATTEMPTS=0          # 0 = unlimited
MAX_SECONDS=0           # 0 = unlimited
SEED=0                  # 0 = get_chm draws a fresh seed each attempt
THREADS=""              # empty = get_chm uses all cores
KEEP_REJECTS=0
WRITE_LAMBDA=0          # also save the invariant set of each keeper
VERIFY=0                # re-check each keeper with summary
QUIET=0
BINDIR=""               # empty = directory of this script

usage() {
    cat <<'EOF'
usage: search_chm.sh [options]

  -k, --count K        matrices to collect                     (default 1)
  -l, --lmin L1        minimum #L, inclusive                   (default 0)
  -u, --lmax L2        maximum #L, inclusive; -1 = no bound     (default -1)
  -N, --size N         order of the matrix                     (default 11)
  -s, --symmetry X     s | h | none                            (default s)
  -m, --method M       sinkhorn | rwcp                         (default sinkhorn)
  -i, --iterations M   passed to get_chm --iterations           (default 1e5)
      --lambda-eps E   passed to lambda --eps
  -o, --outdir DIR     where keepers go                (default ../CHM_NUMERICAL)
      --max-attempts A give up after A attempts, 0 = never      (default 0)
      --max-seconds S  give up after S seconds, 0 = never       (default 0)
      --seed S         base seed; attempt n uses S+n, so a run
                       is reproducible. 0 = fresh OS seed each  (default 0)
  -t, --threads T      passed to get_chm --threads. restarts race
                       across threads, so bit-exact replay of a
                       --seed run needs --threads 1  (default: all cores)
      --keep-rejects   move non-matching matrices to DIR/rejects
                       instead of deleting them
      --write-lambda   also write each keeper's invariant set
      --verify         re-check every keeper with summary
      --bindir DIR     directory holding get_chm/lambda/summary
  -q, --quiet          only print the final summary
  -h, --help           this message

Every run collects into one directory, ../CHM_NUMERICAL, rather than a
separate chm_N<N>_<symmetry> per (size, symmetry).  Each keeper is moved
there, and one line is appended to CHM_NUMERICAL.tsv next to this script --
the matrices may be deleted, the catalog is the lasting record.  Its rows are
in chronological order:

    file  N  symmetry  #L  ud  q  attempt  seed

exit: 0 = K matrices collected
      1 = a limit was reached first (partial results are kept and reported)
      2 = bad usage, or a search that cannot succeed
EOF
}

# ------------------------------------------------------------------- parse ---
need_value() {
    if [ "$2" -eq 0 ]; then
        printf 'search_chm.sh: %s needs a value\n' "$1" >&2
        exit 2
    fi
}

while [ $# -gt 0 ]; do
    case "$1" in
        -k|--count)      need_value "$1" $(($# - 1)); COUNT=$2;        shift 2 ;;
        -l|--lmin)       need_value "$1" $(($# - 1)); LMIN=$2;         shift 2 ;;
        -u|--lmax)       need_value "$1" $(($# - 1)); LMAX=$2;         shift 2 ;;
        -N|--size)       need_value "$1" $(($# - 1)); SIZE=$2;         shift 2 ;;
        -s|--symmetry)   need_value "$1" $(($# - 1)); SYMMETRY=$2;     shift 2 ;;
        -m|--method)     need_value "$1" $(($# - 1)); METHOD=$2;       shift 2 ;;
        -i|--iterations) need_value "$1" $(($# - 1)); ITERATIONS=$2;   shift 2 ;;
        --lambda-eps)    need_value "$1" $(($# - 1)); LAMBDA_EPS=$2;   shift 2 ;;
        -o|--outdir)     need_value "$1" $(($# - 1)); OUTDIR=$2;       shift 2 ;;
        --max-attempts)  need_value "$1" $(($# - 1)); MAX_ATTEMPTS=$2; shift 2 ;;
        --max-seconds)   need_value "$1" $(($# - 1)); MAX_SECONDS=$2;  shift 2 ;;
        --seed)          need_value "$1" $(($# - 1)); SEED=$2;         shift 2 ;;
        -t|--threads)    need_value "$1" $(($# - 1)); THREADS=$2;      shift 2 ;;
        --bindir)        need_value "$1" $(($# - 1)); BINDIR=$2;       shift 2 ;;
        --keep-rejects)  KEEP_REJECTS=1; shift ;;
        --write-lambda)  WRITE_LAMBDA=1; shift ;;
        --verify)        VERIFY=1;       shift ;;
        -q|--quiet)      QUIET=1;        shift ;;
        -h|--help)       usage; exit 0 ;;
        *) printf 'search_chm.sh: unknown option %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

is_uint() { case "$1" in ''|*[!0-9]*) return 1 ;; *) return 0 ;; esac; }

for pair in "COUNT $COUNT" "LMIN $LMIN" "SIZE $SIZE" \
            "MAX_ATTEMPTS $MAX_ATTEMPTS" "MAX_SECONDS $MAX_SECONDS" "SEED $SEED"; do
    set -- $pair
    if ! is_uint "$2"; then
        printf 'search_chm.sh: %s must be a non-negative integer, got "%s"\n' "$1" "$2" >&2
        exit 2
    fi
done
if [ -n "$THREADS" ] && ! is_uint "$THREADS"; then
    printf 'search_chm.sh: --threads must be a non-negative integer\n' >&2
    exit 2
fi
if [ "$LMAX" != "-1" ] && ! is_uint "$LMAX"; then
    printf 'search_chm.sh: --lmax must be a non-negative integer or -1\n' >&2
    exit 2
fi
if [ "$LMAX" != "-1" ] && [ "$LMAX" -lt "$LMIN" ]; then
    printf 'search_chm.sh: --lmax (%s) is below --lmin (%s)\n' "$LMAX" "$LMIN" >&2
    exit 2
fi
if [ "$COUNT" -lt 1 ]; then
    printf 'search_chm.sh: --count must be at least 1\n' >&2
    exit 2
fi
case "$SYMMETRY" in s|h|none) ;; *)
    printf 'search_chm.sh: --symmetry must be s, h or none\n' >&2; exit 2 ;;
esac

# --------------------------------------------------------------- locate ------
if [ -z "$BINDIR" ]; then
    BINDIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
fi
GET_CHM="$BINDIR/get_chm"
LAMBDA="$BINDIR/lambda"
SUMMARY="$BINDIR/summary"
DEFECT="$BINDIR/defect"

for b in "$GET_CHM" "$LAMBDA"; do
    if [ ! -x "$b" ]; then
        printf 'search_chm.sh: %s not found or not executable.\n' "$b" >&2
        printf '  build it first, or pass --bindir; see COMPILE.txt\n' >&2
        exit 2
    fi
done
[ "$VERIFY" -eq 1 ] && [ ! -x "$SUMMARY" ] && { printf 'search_chm.sh: --verify needs %s\n' "$SUMMARY" >&2; exit 2; }

[ -z "$OUTDIR" ] && OUTDIR="$(cd -- "$BINDIR" && pwd)/../CHM_NUMERICAL"
mkdir -p -- "$OUTDIR" || exit 2
[ "$KEEP_REJECTS" -eq 1 ] && mkdir -p -- "$OUTDIR/rejects"
INDEX="$BINDIR/CHM_NUMERICAL.tsv"
if [ ! -s "$INDEX" ]; then
    printf '#file\tN\tsymmetry\t#L\tud\tq\tattempt\tseed\n' > "$INDEX"
fi

# A request that provably has no solution should fail now, not after an hour.
# get_chm's exit codes separate the two cases cleanly: 2 = bad usage or a
# combination with no solution (a Hermitian CHM of odd non-square order, say),
# 1 = the search merely ran out of budget, which is expected of this probe.
probe_err=$(mktemp) || exit 2
"$GET_CHM" --size "$SIZE" --symmetry "$SYMMETRY" --method "$METHOD" \
           --restarts 1 --iterations 1 --threads 1 \
           --out "$probe_err.m" >/dev/null 2>"$probe_err"
probe_rc=$?
rm -f -- "$probe_err.m"
if [ "$probe_rc" -eq 2 ]; then
    cat -- "$probe_err" >&2
    rm -f -- "$probe_err"
    exit 2
fi
rm -f -- "$probe_err"

# ---------------------------------------------------------------- report -----
FOUND=0
ATTEMPT=0
START=$SECONDS
declare -a KEPT=()

say() { [ "$QUIET" -eq 1 ] || printf '%s\n' "$*"; }

report() {
    local el=$((SECONDS - START))
    printf '\n'
    printf '%s\n' "-----------------------------------------------------------"
    printf 'collected %d of %d requested matrices\n' "$FOUND" "$COUNT"
    printf 'N = %s, symmetry = %s, method = %s, #L in [%s, %s]\n' \
           "$SIZE" "$SYMMETRY" "$METHOD" "$LMIN" \
           "$([ "$LMAX" = "-1" ] && echo 'inf' || echo "$LMAX")"
    printf 'attempts = %d, elapsed = %ds' "$ATTEMPT" "$el"
    if [ "$ATTEMPT" -gt 0 ]; then
        printf ', hit rate = %d%%' $((100 * FOUND / ATTEMPT))
    fi
    printf '\n'
    if [ "${#KEPT[@]}" -gt 0 ]; then
        printf 'kept in %s:\n' "$OUTDIR"
        for f in "${KEPT[@]}"; do printf '  %s\n' "$f"; done
    fi
    printf 'index: %s\n' "$INDEX"
    printf '%s\n' "-----------------------------------------------------------"
}

interrupted=0
on_int() { interrupted=1; }
trap on_int INT TERM

# ----------------------------------------------------------------- search ----
say "searching for $COUNT matrix(es): N = $SIZE, symmetry = $SYMMETRY, method = $METHOD"
say "keeping #L in [$LMIN, $([ "$LMAX" = "-1" ] && echo 'inf' || echo "$LMAX")]  ->  $OUTDIR"
say ""

while [ "$FOUND" -lt "$COUNT" ]; do
    [ "$interrupted" -eq 1 ] && { say ""; say "interrupted"; break; }
    if [ "$MAX_ATTEMPTS" -gt 0 ] && [ "$ATTEMPT" -ge "$MAX_ATTEMPTS" ]; then
        say ""; say "attempt limit reached ($MAX_ATTEMPTS)"; break
    fi
    if [ "$MAX_SECONDS" -gt 0 ] && [ $((SECONDS - START)) -ge "$MAX_SECONDS" ]; then
        say ""; say "time limit reached (${MAX_SECONDS}s)"; break
    fi

    ATTEMPT=$((ATTEMPT + 1))

    # attempt seed: 0 keeps get_chm's own OS seed, otherwise SEED + attempt so
    # that the whole run can be replayed exactly
    if [ "$SEED" -eq 0 ]; then seed_arg=(); this_seed="os"
    else seed_arg=(--seed $((SEED + ATTEMPT))); this_seed=$((SEED + ATTEMPT)); fi
    if [ -n "$THREADS" ]; then thread_arg=(--threads "$THREADS"); else thread_arg=(); fi

    FILE=$("$GET_CHM" --size "$SIZE" --symmetry "$SYMMETRY" --method "$METHOD" \
                      --iterations "$ITERATIONS" "${seed_arg[@]}" \
                      "${thread_arg[@]}" --quiet 2>/dev/null)
    rc=$?
    if [ $rc -ne 0 ] || [ -z "$FILE" ] || [ ! -f "$FILE" ]; then
        say "$(printf '[%3d/%3d] attempt %-5d get_chm failed (rc=%d)' \
              "$FOUND" "$COUNT" "$ATTEMPT" "$rc")"
        [ $rc -eq 2 ] && { say "  usage error -- aborting"; report; exit 2; }
        continue
    fi

    # ---- the filter --------------------------------------------------------
    if [ -n "$LAMBDA_EPS" ]; then
        L=$("$LAMBDA" --file "$FILE" --eps "$LAMBDA_EPS" --bare --no-write 2>/dev/null)
    else
        L=$("$LAMBDA" --file "$FILE" --bare --no-write 2>/dev/null)
    fi
    if ! is_uint "${L:-}"; then
        say "$(printf '[%3d/%3d] attempt %-5d lambda failed on %s' \
              "$FOUND" "$COUNT" "$ATTEMPT" "$FILE")"
        rm -f -- "$FILE"
        continue
    fi

    keep=1
    [ "$L" -lt "$LMIN" ] && keep=0
    [ "$LMAX" != "-1" ] && [ "$L" -gt "$LMAX" ] && keep=0

    if [ "$keep" -eq 0 ]; then
        say "$(printf '[%3d/%3d] attempt %-5d #L = %-8s reject' \
              "$FOUND" "$COUNT" "$ATTEMPT" "$L")"
        if [ "$KEEP_REJECTS" -eq 1 ]; then mv -- "$FILE" "$OUTDIR/rejects/"
        else rm -f -- "$FILE"; fi
        continue
    fi

    # ---- a keeper ----------------------------------------------------------
    if [ "$VERIFY" -eq 1 ] && ! "$SUMMARY" --file "$FILE" >/dev/null 2>&1; then
        say "$(printf '[%3d/%3d] attempt %-5d #L = %-8s FAILED VERIFY -- dropped' \
              "$FOUND" "$COUNT" "$ATTEMPT" "$L")"
        rm -f -- "$FILE"
        continue
    fi

    UD="-"
    [ -x "$DEFECT" ] && UD=$("$DEFECT" --file "$FILE" --bare 2>/dev/null || echo '-')
    Q="-"
    if [ -x "$SUMMARY" ]; then
        Q=$("$SUMMARY" --file "$FILE" --butson 2>/dev/null \
            | sed -n 's/^\* this is a BH([0-9]*, \([0-9]*\))$/\1/p')
        [ -z "$Q" ] && Q="-"
    fi

    mv -- "$FILE" "$OUTDIR/" || { say "  cannot move $FILE"; continue; }
    BASE=$(basename -- "$FILE")
    FOUND=$((FOUND + 1))
    KEPT+=("$BASE")

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n' \
           "$BASE" "$SIZE" "$SYMMETRY" "$L" "$UD" "$Q" "$ATTEMPT" "$this_seed" >> "$INDEX"

    if [ "$WRITE_LAMBDA" -eq 1 ]; then
        LOUT="$OUTDIR/lambda_${BASE%.m}.data"
        if [ -n "$LAMBDA_EPS" ]; then
            "$LAMBDA" --file "$OUTDIR/$BASE" --eps "$LAMBDA_EPS" --out "$LOUT" >/dev/null 2>&1
        else
            "$LAMBDA" --file "$OUTDIR/$BASE" --out "$LOUT" >/dev/null 2>&1
        fi
    fi

    say "$(printf '[%3d/%3d] attempt %-5d #L = %-8s ud = %-4s KEEP  %s' \
          "$FOUND" "$COUNT" "$ATTEMPT" "$L" "$UD" "$BASE")"
done

report
[ "$FOUND" -ge "$COUNT" ] && exit 0
exit 1
