# Complex Hadamard matrices — software specification

*Last update: 2026-09-13.*

This document describes the toolbox as it stands; it is not a changelog.

## Contents

- [Conventions](#conventions)
- [Build](#build)
- [Quick start](#quick-start)
- [`get_chm.cpp`](#get_chmcpp) — return a complex Hadamard matrix of given size
- [`defect.cpp`](#defectcpp) — the defect of a unitary matrix
- [`lambda.cpp`](#lambdacpp) — cardinality of the Haagerup invariants
- [`summary.cpp`](#summarycpp) — check whether a matrix is a CHM
- [`check_me.cpp`](#check_mecpp) — monomial equivalence of two matrices
- [`symmetrize.cpp`](#symmetrizecpp) — put a matrix into symmetric form
- [Is every class symmetrisable?](#is-every-class-symmetrisable)
- [Input files accepted](#input-files-accepted)
- [Tolerances](#tolerances)
- [Existence constraints](#existence-constraints)
- [Validation](#validation)
- [Open questions](#open-questions)

---

## Conventions

A complex Hadamard matrix (CHM) of order `N` satisfies

```
|H_jk| = 1              for all j, k
H * H' = N * I_N
```

measured by the two norms used everywhere below:

```
nf = || H*H' - N*I ||_F         (called nh in the Octave scripts)
n1 = ||  |H| - 1    ||_F
```

**Dephased form** — first row and first column are all ones:

```
H = [ 1 1 1 ... 1     ]
    [ 1               ]
    [ 1   CORE        ]
    [ :   (N-1)x(N-1) ]
    [ 1               ]
```

**Core phases** — the real `(N-1)x(N-1)` array `A` with `CORE = exp(2j*pi*A)`.
Phases are normalised to the unit interval, `A(j,k)` in `[0, 1)`. The full matrix
is recovered by

```matlab
>> H = [ones(1, columns(A) + 1); [ones(rows(A), 1), exp(2j*pi*A)]];
```

This dephasing convention, `K = diag(1./H(:,1)) * H * diag(1./H(1,:))`, preserves
both `H = H^T` and `H = H'`, so a symmetry requested on the command line survives
into the saved file.

---

## Build

Dependencies, Ubuntu or Debian:

```sh
sudo apt install build-essential liblapacke-dev libopenblas-dev
```

`liblapacke-dev` provides `lapacke.h` and `-llapacke` (the SVD, `dgesdd` and
`zgesdd`), `libopenblas-dev` a fast BLAS that dispatches on the CPU at run time,
and `libgomp` — OpenMP — ships with g++. Nothing else is needed: no CMake, no
Boost, no Eigen.

Tuned for the target machine — AMD Ryzen 5 5600H, Zen 3, 6 cores / 12 threads,
AVX2 + FMA, no AVX-512, L3 16 MiB:

```sh
for p in get_chm defect lambda summary check_me symmetrize; do
    g++ -std=c++17 -O3 -march=znver3 -mtune=znver3 -fopenmp -flto \
        -funroll-loops -fno-math-errno \
        -o $p $p.cpp chm.cpp -llapacke -llapack -lopenblas -lm
done
```

Portable on any x86-64 Ubuntu machine (`-march=znver3` would SIGILL elsewhere):

```sh
g++ -std=c++17 -O3 -march=x86-64-v3 -fopenmp -flto \
    -o get_chm get_chm.cpp chm.cpp -llapacke -llapack -lopenblas -lm
```

The choices worth knowing:

| `-march=` | what it assumes |
|---|---|
| `x86-64` | the safe baseline, runs on anything x86-64 |
| `x86-64-v2` | SSE4.2 + POPCNT, anything from about 2009 |
| `x86-64-v3` | AVX2 + FMA, Haswell / Zen 1 and later — a good default |
| `znver3` | the Ryzen 5 5600H specifically |
| `native` | whatever CPU is compiling — simplest, if the binary never travels |

Drop `-march`/`-mtune` entirely on a non-x86 host; the sources contain no
intrinsics and no assembly. Little speed is lost either way, because the
expensive step — the SVD — runs inside OpenBLAS, which dispatches on the CPU at
run time whatever this code was compiled for.

`-DNDEBUG` removes the assertions, though none of them sit in a hot path, and
`-fsanitize=address,undefined -O1 -g` gives a debug build. Dropping `-fopenmp`
still compiles and runs: the restart loop and the Haagerup scan simply become
single-threaded.

### `-march=znver3` is not portable and not a mere tuning hint

It enables instruction sets only AMD Zen 3 has — `sse4a`, `sha`, `clwb`,
`clzero`, `mwaitx`, `rdpid`, `vaes`, `vpclmulqdq` — and g++ never looks at the CPU
it is running on, so the build succeeds anywhere and the binary may then die with
SIGILL. On an Intel i7-7700 (Kaby Lake), which has `avx2`, `fma`, `bmi2` and
`f16c` but none of that list, use **`-march=native`** or **`-march=x86-64-v3`**;
the latter satisfies both machines. To see what a given choice actually permits:

```sh
g++ -march=znver3    -Q --help=target | grep -E 'sse4a|sha|clwb|vaes'
g++ -march=x86-64-v3 -Q --help=target | grep -E 'sse4a|sha|clwb|vaes'
```

The first prints `[enabled]` for all four, the second `[disabled]` — while
keeping `-mavx2`, `-mfma` and `-mbmi2` on, which is exactly what the i7-7700 has.

### Extract into an empty directory

`chm.hpp`, `chm.cpp` and the programs are one versioned unit and code
moves between them between revisions, so unpacking a new archive over an old tree
leaves a mixture that does not compile. Each `.cpp` declares `CHM_PROGRAM_ABI`
before including `chm.hpp`, which declares `CHM_ABI`, and a mismatch in either
direction stops the build with one line naming the cause instead of a page of
overload-resolution notes:

```
#error "MIXED SOURCE TREE: ... extract one archive into it"
```

The case this guards against is real: before 2026-09-08 `get_chm.cpp` carried its
own file-static `exe_dir()` and `run_int()`, which `chm.hpp` now declares, and
with `using namespace chm;` in force every call to them became ambiguous.

### Threads

Restarts of the search run in parallel, and OpenBLAS is pinned to one thread
inside that region so the two thread pools do not oversubscribe the cores. The
default is one thread per hardware thread — 12 on the Ryzen 5 5600H, 8 on a
4-core i7-7700 — and nothing needs changing between the two:

```sh
./get_chm --size 12 --threads 6
OMP_NUM_THREADS=6 ./get_chm --size 12
```

### Files

| File | Role |
|---|---|
| `chm.hpp` | shared interface (header only) |
| `chm.cpp` | shared library (compiled into every program) |
| `cli.hpp` | command-line helper (header only) |
| `get_chm.cpp` `defect.cpp` `lambda.cpp` `summary.cpp` | the original four programs |
| `check_me.cpp` `symmetrize.cpp` | monomial equivalence, and symmetric form |

**Build `defect` and `lambda` before using `get_chm`** — `get_chm` calls both of
them to name its output file (see [Output](#output)). They are looked for next to
the `get_chm` binary itself, not in the current directory, so the toolbox works
from anywhere; if either is missing the neutral fallback name is used and nothing
else changes.

---

## Quick start

```console
$ ./get_chm                                # CHM(7), sinkhorn, defaults
$ F=$(./get_chm --size 9 --quiet)          # e.g. SH_9_0_681_20260904T224151.m
$ G=$(./get_chm --size 9 --quiet)
$ ./summary --file "$F" --all
$ ./defect  --file "$F"
$ ./lambda  --file "$F"
$ ./check_me "$F" "$G"                     # 0 = equivalent, 1 = not
$ ./symmetrize --file "$F"                 # 0 = a symmetric form exists
$ ./check_me "$F" "$F" --count             # |Aut(H)/U(1)|
$ ./get_chm --help                         # every flag, with its default
```

The name `get_chm` writes is `xH_N_d_L_<timestamp>.m`: `x` is `S`, `H` or `Y`
for symmetric, Hermitian or neither, `N` the order, then the dephased defect and
the number of Haagerup invariants, obtained by calling `./defect` and
`./lambda`. If either helper cannot be run the name falls back to
`A_<timestamp>.m`, and an explicit `--out FILE` is always used verbatim.
`symmetrize` names its result the same way, reusing `d` and `L` from the input
file name when it carries them — both are invariants of the equivalence class.

---

## `get_chm.cpp`

*Return a complex Hadamard matrix of given size.*

### Flags

| Flag | Meaning | Default |
|---|---|---|
| `--size N` | order of the matrix | `7` |
| `--method M` | `sinkhorn` = Sinkhorn alternating-projection algorithm; `rwcp` = random walk over core phases | `sinkhorn` |
| `--iterations M` | maximum iterations per restart; both `1e+4` and `10000` are accepted. `--iteration` is kept as an alias for the singular spelling | `1e+4` |
| `--symmetry X` | `s` = symmetric (`H = H^T`), `h` = Hermitian (`H = H'`), `none` | `none` |
| `--seed S` | PRNG seed; `0` draws one from the OS | `0` |
| `--restarts R` | maximum restarts, `0` = until success | `0` |
| `--eps E` | convergence threshold on `nf` and `n1` | `1e-13` |
| `--digits D` | significant digits written to the file | `17` |
| `--threads T` | parallel restart threads, `0` = all cores | `0` |
| `--kicks K` | phase perturbations tried before abandoning a restart, `0` = restart at once | `8` |
| `--max-seconds S` | wall-clock budget, `0` = unlimited | `0` |
| `--no-polish` | skip the exact-coordinate-descent stage that sharpens a near-miss sinkhorn iterate | — |
| `--out FILE` | output path, used verbatim; without it the name is built from the matrix itself, `xH_N_d_L_YYYYMMDDThhmmss.m` | — |
| `--quiet` | print the output file name only | — |
| `--verbose` | one dot per failed restart, on stderr | — |
| `--help`, `-h` | every flag with its default | — |

Default setting, returns a CHM of size 7 using the default method:

```console
$ ./get_chm
```

Other examples:

```console
$ ./get_chm --size N
$ ./get_chm --size N --method sinkhorn
$ ./get_chm --size N --method rwcp
$ ./get_chm --size 6 --symmetry s
$ ./get_chm --size 8 --symmetry h
$ ./get_chm --size 16 --iterations 3e+5 --seed 20260904
```

### Exit codes

| Code | Meaning |
|---|---|
| 0 | a CHM was found and written |
| 1 | the search failed within the given budget |
| 2 | bad usage, or a request that provably has no solution |

### Output

The matrix is written in the dephased form, as the `(N-1)x(N-1)` array of core
phases. The file is named

```
xH_N_d_L_YYYYMMDDThhmmss.m
```

| Field | Meaning |
|---|---|
| `x` | `S` symmetric, `\|\| H - H^T \|\|_F <= 1e-7`; `H` Hermitian, `\|\| H - H' \|\|_F <= 1e-7`; `Y` neither |
| `N` | the order of the matrix |
| `d` | the dephased defect, obtained by running `./defect --bare` |
| `L` | the number of Haagerup invariants, from `./lambda --bare --no-write` |

So `SH_11_0_139_20260904T162711.m` is a symmetric CHM of order 11, defect 0,
lying in the `#Lambda = 139` stratum. The fields run from the coarsest to the
finest, so a glob is already a query and no bookkeeping is needed on the side:

```sh
ls SH_11_0_*        # every isolated symmetric CHM(11) collected so far
ls ?H_9_*           # everything of order 9, whatever its symmetry
ls SH_11_0_139_*    # one stratum
```

`N` is redundant — it is (size of the array) + 1 — but it is what one filters on
most often, and reading it off the name beats opening the file.

The fields sort lexicographically, not numerically, so `SH_11_` precedes `SH_9_`
in an `ls` listing. `ls -v`, or a glob per order, gives the natural order.

The two helpers are separate programs and each takes a few seconds at moderate
`N`, so the run reports what it is waiting for:

```console
$ ./get_chm --size 9 --symmetry s
* now d (the defect) is being calculated ...
*   d = 0
* now L (the Haagerup invariants) is being calculated ...
*   L = 681
* CHM(9) found: method = sinkhorn, symmetry = s
  ...
* 8 x 8 array of core phases written to SH_9_0_681_20260904T224151.m (17 digits)
```

The matrix is saved first, under the neutral name, so a crash or an interrupt in
the middle of the two calls cannot lose it; the final file is then written afresh
from the same matrix in memory and the neutral one is removed. It is written
rather than renamed because Octave resolves a function by the file's base name
and the `function A = ...` line has to follow the file — a renamed file would
carry the old identifier and could not be called.

**No sanity check is made beyond "did the helper run and print one
integer".** If `./defect` or `./lambda` is missing, is not executable or fails for
any reason, the file simply keeps the neutral name

```
A_YYYYMMDDThhmmss.m
```

and the run says so:

```
*   could not run ./lambda
*   keeping the neutral name A_20260904T223736.m
```

An explicit `--out FILE` is used exactly as given and switches the whole naming
step off; `--quiet` keeps stdout to the one line holding the final file name and
sends the progress lines to stderr, so

```sh
F=$(./get_chm --size 9 --symmetry s --quiet)
```

still works and `$F` is the renamed file. As before, a collision within the same
second appends `_1`, `_2`, … to the name.

Example. `./get_chm` at the default `N = 7` saves a 6×6 array of real numbers.
This is the Octave m-file:

```matlab
function A = YH_7_4_451_20260904T223739
% ------------------------------------------------------------------------------
% 20260904T223739
% core phases
% usage: >> A = YH_7_4_451_20260904T223739; H = [ones(1, columns(A) + 1); [ones(rows(A), 1), exp(2j*pi*A)]];
% phases are normalised to the unit interval: A(j,k) in [0, 1)
% N = 7, method = sinkhorn, symmetry = none                     <- new
% seed = 7056868887342391694, iterations = 156, restarts = 2     <- new
% nh = || H*H' - N*I ||_F = 9.496e-13                            <- new
% n1 = ||  |H| - 1     ||_F = 2.022e-13                          <- new
% ------------------------------------------------------------------------------
    A = [
        0.49823792834 0.75119220948 ... ;
        0.13094871129 0.90226002550 ... ;
        ...
    ];
end
```

**The function name equals the file's base name**, because Octave and
MATLAB resolve a function by file name: a file `YH_7_4_451_20260904T223739.m`
declaring `function A = core_phases` cannot be called at all. This is also why
the name starts with a letter and not with the defect: an Octave identifier may
not begin with a digit, and `A_...`, `SH_...`, `HH_...`, `YH_...` all do.

**Rows of the array are separated by a semicolon**; a single flat list of
`(N-1)^2` numbers would parse as a row vector, not as a square array.

Note that the phases are normalised to the unit interval.

### Practical iteration counts

*(measured, 12 threads)*

`sinkhorn` is an alternating projection: it converges linearly and its composite
map has fixed points that are not CHM. The default `M = 1e+4` is ample at the
default `N = 7` but not beyond about `N = 12`, where the iterate parks at
`n1 ~ 1e-2` until a phase kick or a restart moves it to another basin.

| `N` | `--iterations` | typical restarts |
|---|---|---|
| 5 … 10 | `1e+4` | 1 … 10 |
| 12 | `1e+5` | tens |
| 16 | `3e+5` | a few |
| 20 | `1e+6` | a few |
| 24 and up | `1e+6` and up, or `--method rwcp` for small `N` only | |

When the search fails the program reports the closest approach reached and, if it
came close, suggests raising `--iterations`. `--max-seconds` keeps an unattended
run from spinning for ever.

`rwcp` walks over the `(N-1)^2` core phases directly, so `|H_jk| = 1` and the
dephased form hold exactly at every step and only `nf` has to be driven to zero:

```
F(A) = nf^2 = sum_{a != b} |G_ab|^2 ,   G = H*H'
```

(the diagonal of `G` is exactly `N` because every entry is unimodular). Moving a
single phase changes one entry of `H`, hence one row and one column of `G`, so
`F` updates in `O(N)` instead of `O(N^3)`. Each step draws a random core
coordinate and either perturbs it by an adaptive step or jumps to its exact
minimiser — with `u = H_rc` on the unit circle the `r`-part of `F` is
`const + 4*Re(u*Z)`, minimal at `u = -conj(Z)/|Z|`, also `O(N)`. Proposals that
would raise `F` are rejected; a stalled walk is kicked by re-randomising a few
phases. `rwcp` is fast and exact for small `N` but its landscape traps it past
about `N = 10`; `sinkhorn` is the better choice for larger orders.

---

## `defect.cpp`

*Calculate the defect of a unitary matrix.*

The dephased defect

```
d(U) = (N-1)^2 - rank(R)
```

where `R` is the real `2*tau x N^2` system matrix (`tau = N(N-1)/2`) of the
linearised phase-perturbation conditions, split into real and imaginary parts.
See W. Tadej, K. Życzkowski, *Linear Algebra Appl.* **429**, 447–481 (2008).

### Flags

| Flag | Meaning | Default |
|---|---|---|
| `--file FILE` | phase array written by `get_chm`; **mandatory** | — |
| `--layout L` | `core` = FILE holds the `(N-1)x(N-1)` core array; `full` = FILE holds the full `NxN` phase array; `auto` = whichever is closer to a CHM | `core` |
| `--method M` | `R` = exact rank of `R`; `S` = rank as `#{ singular values of R > --tol }`; `T` = tangent-space method | `S` |
| `--tol T` | singular-value threshold for method `S` | `1e-8` |
| `--bare` | print the integer only | — |

```console
$ ./defect --file SH_9_0_681_20260904T224151.m
$ ./defect --file SH_9_0_681_20260904T224151.m --method S --tol 1e-10
$ ./defect --file SH_9_0_681_20260904T224151.m --bare
```

Method `S` is the default because a numerically obtained matrix is never exact.
An entry of `1.000000` arriving as `0.999998` leaves the plain rank
untrustworthy, while a singular-value threshold is unmoved by it.

**Methods `R` and `T` always agree.** The tangent-space matrix `T`, with
entries `T(k-th block, pair) = Re(H_k .* conj(H_l))` and `Im(...)`, is exactly
the transpose of `R` entry by entry, and a matrix and its transpose share their
rank. One code path serves both.

---

## `lambda.cpp`

*Calculate the cardinality of the (unique) Haagerup invariants.*

```
Lambda(H) = { H_ij * H_kl * conj(H_il) * conj(H_kj) : i,j,k,l = 1..N }
```

Additionally it writes the unique set of invariants to
`lambda_YYYYMMDDThhmmss.data`. Each line contains two real numbers, the real and
the imaginary part, separated by a tab; lines are sorted by argument. Example:

```
0.283984240234098	0.987498232040958
0.483183318451732	0.526348767823444
1.000000000000043	0.000000000003312
...
```

### Flags

| Flag | Meaning | Default |
|---|---|---|
| `--file FILE` | phase array written by `get_chm`; **mandatory** | — |
| `--layout L` | `core` \| `full` \| `auto` | `core` |
| `--eps E` | two invariants count as equal when their distance is below `E` | `1e-8` |
| `--out FILE` | output path | `lambda_YYYYMMDDThhmmss.data` |
| `--no-write` | report the cardinality without writing the set | — |
| `--bare` | print the cardinality only | — |

```console
$ ./lambda --file SH_9_0_681_20260904T224151.m
$ ./lambda --file SH_9_0_681_20260904T224151.m --eps 1e-10 --bare
```

**Implementation note.** With `r_j = H_ij / H_kj` the quadruple invariant
factorises as `r_j / r_l`, so each row pair `(i,k)` needs `N` ratios and not
`N^2` four-way products; pairs `(i,k)` and `(k,i)` give reciprocal sets, and
`(j,l)` already runs over both orders, so only `i < k` is visited. Uniqueness is
decided by a hash grid of cell size `eps` (9 cells probed per candidate), which
is `O(M)` expected. A pairwise scan would be `O(M^2)`, and with `M = O(N^4)`
candidates that is unusable past about `N = 12`.
The scan is OpenMP-parallel over row pairs.

---

## `summary.cpp`

*Check whether a given matrix is a CHM.*

### Flags

| Flag | Meaning | Default |
|---|---|---|
| `--file FILE` | the matrix; **mandatory** | — |
| `--symmetry` | check whether the matrix is symmetric or Hermitian | — |
| `--defect` | calculate the defect | — |
| `--butson` | check whether the matrix belongs to some Butson class | — |
| `--lambda` | calculate the cardinality of the Haagerup invariants | — |
| `--SL3` | calculate the linear entropies | — |
| `--all` | everything above | — |
| `--layout L` | `core` \| `full` \| `auto` | `core` |
| `--eps E` | CHM threshold on `nf` and `n1` | `1e-7` |
| `--tol T` | singular-value threshold for the defect | `1e-8` |
| `--lambda-eps E` | uniqueness threshold for the invariants | `1e-8` |
| `--qmax Q` | largest Butson exponent tested | `10000` |
| `--butson-tol T` | tolerance of the Butson test | `1e-8` |
| `--log-form` | print the LOG-form matrix when `H` is Butson | — |

```console
$ ./summary
no arguments prints warning

$ ./summary --file xH_N_d_L_YYYYMMDDThhmmss.m
$ ./summary --file xH_N_d_L_YYYYMMDDThhmmss.m --symmetry
```

Checking symmetry is a default option; other options must be written explicitly,
in any order.

```console
$ ./summary --file xH_N_d_L_YYYYMMDDThhmmss.m --defect
$ ./summary --file xH_N_d_L_YYYYMMDDThhmmss.m --defect --SL3
$ ./summary --file xH_N_d_L_YYYYMMDDThhmmss.m --lambda --defect --symmetry
$ ./summary --file xH_N_d_L_YYYYMMDDThhmmss.m --all --log-form
```

It takes the `.m` file produced by `get_chm.cpp` and checks:

- is this a CHM (`H = exp(2j*pi*A)`; `nf = || H*H' - N*I ||_F`; `n1 = || |H| - 1 ||_F`)
- is this a Butson-type matrix
- is this a symmetric or Hermitian matrix
- the cardinality of the Haagerup invariants
- the value of the defect
- the values of the linear entropies of `H`, `H^R` and `H^G`

If `nf > 1e-7` or `n1 > 1e-7` it prints `* this is not a CHM!` and ends with
**exit code 1**. Otherwise it prints:

```
* this is a CHM (nf = $nf, n1 = $n1)
* this is a BH($N, $q) [or this is probably not a BH-matrix for any q]
* this is a symmetric matrix [or nothing if it is not symmetric]
* this is a Hermitian matrix [or nothing if it is not Hermitian]
* #L = $L
* ud = $d [unitary defect]
* SL(H) = $s1 $s2 $s3 [where $sj is the entropy of H, H^R or H^G]
```

Variables `$...` are set appropriately.

**Only the CHM test and the symmetry lines are unconditional**; the BH,
`#L`, `ud` and SL lines appear when their flag (or `--all`) is given. This
resolves the earlier conflict between the flag list, which made `--butson`
explicit, and the output list, which showed the BH line always.

**`SL(H)` needs `N = d*d`.** The reshuffling `H^R` and the partial
transpose `H^G` are defined only for a bipartite `d x d` split, so for any other
order the line reads

```
* SL(H) = n/a  (H^R and H^G need N = d*d; here N = 7)
```

`H^G` is the partial transpose over the second subsystem. Transposing the first
instead differs by a global transpose, which leaves every singular value and
therefore SL unchanged, so either convention gives the same triplet.
After the triplet, `R-dual matrix!`, `G-dual matrix!` or `2-unitary matrix!` is
printed when the corresponding entropies equal 1.

### Exit codes

| Code | Meaning |
|---|---|
| 0 | the matrix is a CHM |
| 1 | the matrix is not a CHM |
| 2 | bad usage, or the file cannot be read |

---

## `check_me.cpp`

*(`me` = monomial equivalence)*

*Decide whether two complex Hadamard matrices are monomially equivalent.*

Two CHM are monomially equivalent when

```
H2 = D1 P1 H1 P2 D2
```

with `D1`, `D2` diagonal unitary and `P1`, `P2` permutations. This is the
equivalence under which CHM are classified, and `#Lambda` alone does not decide
it: matrices sharing a stratum are routinely inequivalent, so the question has to
be settled matrix by matrix.

**How it decides.** Dephasing `H` at one of its entries `(a,b)`,

```
dephase_{a,b}(H)_{kl} = H_kl * H_ab / (H_al * H_kb)
```

turns row `a` and column `b` into ones and thereby fixes `D1` and `D2`
completely: the diagonal factors are no longer free, only the two permutations
are. Every dephased representative of the class of `H2` is therefore

```
P_row . dephase_{a,b}(H2) . P_col
```

for one of the `N*N` choices of base `(a,b)`. The program walks all `N*N` bases
and, for each, searches for the permutations explicitly: rows are assigned one at
a time, and a branch is cut as soon as the sorted multiset of column labels of
the partial matrix can no longer be matched. Phases are first replaced by
codebook labels, shared between both matrices, so the comparisons are on
integers.

The search is exhaustive, so a negative answer is a proof up to the phase
tolerance, and a positive answer comes with the witnesses `sigma` and `tau`,
which are then verified on the actual phases — the printed mismatch is that
check.

### Usage

```console
$ ./check_me FILE1 FILE2 [--tol T] [--quiet]
```

Both files are the `.m` core-phase files written by `get_chm`, taken directly,
with no conversion step:

```console
$ ./check_me SH_11_0_139_20260904T162711.m SH_11_0_139_20260905T081530.m
```

Older two-column `.hp` files are still accepted, and so is any file the other
programs read; see [Input files accepted](#input-files-accepted).

| Flag | Meaning | Default |
|---|---|---|
| `--tol T` | two phases count as equal below `T` | `1e-9` |
| `--count` | count every equivalence instead of stopping at the first; with the same file twice this is `\|Aut(H)/U(1)\|` | — |
| `--quiet` | print `EQUIVALENT` / `INEQUIVALENT` and nothing else | — |
| `--help`, `-h` | usage and defaults | — |

**`--count` and the automorphism group.** The columns of a CHM are pairwise
non-proportional, so for a given base and row assignment the column permutation
is unique if it exists at all — one equivalence per hit, and the total is exact.
Run on one matrix against itself it returns the order of the monomial
automorphism group modulo the scalars, which is the group that decides the
symmetrisation question (see [`symmetrize.cpp`](#symmetrizecpp)). Checked against
four known values:

```console
$ ./check_me F_p.m F_p.m --count
```

| Matrix | `--count` | Known value |
|---|---|---|
| `F_p`, `p` = 5, 7, 11 | 100, 294, 1210 | `p^2 (p-1)` |
| the Hadamard matrix of order 12 | 95040 | `\|M_12\|`, the Mathieu group |
| the Hadamard matrix of order 8 | 10752 | 10752 |
| the Sylvester matrix of order 16 | 5160960 | 5160960 |

### Exit codes

| Code | Meaning |
|---|---|
| 0 | the two matrices are monomially equivalent |
| 1 | they are not |
| 2 | bad usage, or a file cannot be read |

```console
$ ./check_me A.m B.m > /dev/null && echo same || echo different
```

**On the tolerance.** `get_chm` writes 17 digits of a matrix converged to
`nf ~ 1e-13`, so its phases are meaningful to about `1e-13` and `1e-9` is the
right default — comfortably above the noise, far below the gap between distinct
phases (reported as "smallest gap" on every run, typically `1e-4 .. 1e-6`).
Asking for much less is a mistake and not a refinement: comparing
double-precision input at `1e-25` reports inequivalence for a matrix and its own
transpose. Only when a genuinely refined matrix is at hand — one whose residual
has been driven far below double precision, by whatever means — is a smaller
`--tol` meaningful. If the smallest gap ever approaches the tolerance the run says so,
because the codebook is then no longer trustworthy.

**Cost.** The base loop is `N^2` and each base is a backtracking search over
row assignments; the column-multiset cut makes it finish in tens of nodes rather
than `N!` in every case tried up to `N = 11`. A positive answer usually stops at
the first base.

### Example output

```console
$ ./check_me circ_full_L11.m circ_full_L11_T.m
* circ_full_L11.m     N = 11, nh = 1.010e-14, n1 = 3.511e-16
* circ_full_L11_T.m   N = 11, nh = 1.010e-14, n1 = 3.511e-16
* codebook: 11 distinct phases, smallest gap 9.091e-02 (tolerance 1e-09)
* base of circ_full_L11_T.m: (row 0, col 0)
* row permutation:    0 1 2 3 4 5 6 7 8 9
* column permutation: 0 1 2 3 4 5 6 7 8 9
* max phase mismatch after applying them: 6.245e-16
* (11 search nodes)
* EQUIVALENT

$ ./check_me NEW_L139.m NEW_L161.m
* codebook: 299 distinct phases, smallest gap 5.579e-05 (tolerance 1e-09)
* all 121 bases exhausted, 121 search nodes
* INEQUIVALENT
```

Matrices of different order are reported inequivalent, not as an error.

---

## `symmetrize.cpp`

*Put a matrix into symmetric form, if its equivalence class admits one.*

```console
$ ./symmetrize --file xH_N_d_L_YYYYMMDDThhmmss.m
```

Decides whether the monomial equivalence class of `H` contains a representative
`S` with `S = S^T`, and writes one out when it does. The search is exhaustive, so
a negative answer is a proof up to the phase tolerance and not a failure to find.

### The reduction that makes it cheap

Dephasing at an entry `(a,b)`,

```
K = dephase_{a,b}(H) ,     K_kl = H_kl * H_ab / (H_al * H_kb)
```

cancels both diagonal factors: for `H' = D1 P H Q D2` one gets, entry by entry,

```
dephase_{a,b}(H') = P . dephase_{p(a),q(b)}(H) . Q                      (*)
```

with the `d`'s gone. So every dephased member of the class is a row and column
permutation of one of the `N*N` matrices `dephase_{a,b}(H)`.

Now let `S = S^T` lie in the class. Dephasing `S` at a **diagonal** entry keeps it
symmetric, because `S_al = S_la` and `S_ka = S_ak` give

```
dephase_{a,a}(S)_kl = S_kl * S_aa / (S_al * S_ka) = dephase_{a,a}(S)_lk
```

so by `(*)` that symmetric matrix is `P K Q` for some base and some `P, Q`.
Writing `P K Q` as `(i,j) -> K_{alpha(i), beta(j)}` and putting
`rho = alpha^-1 beta`, the symmetry condition reads

```
K_{alpha(i), alpha(rho(j))} = K_{alpha(j), alpha(rho(i))}
```

and the left-hand matrix is a **simultaneous** row-and-column relabelling of the
matrix `K'_uv = K_{u, (alpha rho alpha^-1)(v)}` — a relabelling that preserves
symmetry, while `alpha rho alpha^-1` runs over all permutations as `rho` does.
Hence

```
the class of H contains a symmetric matrix
    <=>  for some base (a,b) and some permutation pi,
         K_{i, pi(j)} = K_{j, pi(i)}  for all i, j,   K = dephase_{a,b}(H)
```

One permutation per base, not two, and the witness is `S = K P_pi`, already
symmetric and already a CHM. Two facts prune what is left:

- **`pi(a) = b` is forced.** Row `a` of `K P_pi` is all ones, so by symmetry
  column `a` must be too; column `a` of `K P_pi` is column `pi(a)` of `K`, and
  column `b` is the only all-ones column of `K` — a second one would make two
  columns of `H` proportional, which orthogonality forbids. So the symmetric
  representative is automatically dephased at `(a,a)`.
- **`multiset(row i of K) = multiset(column pi(i) of K)`**, since row `i` of
  `K P_pi` is a permutation of row `i` of `K` and must equal column `i` of
  `K P_pi`. Matching multisets is computed once per base and usually leaves one
  or two candidates per row.

The rest is backtracking: fixing `pi(j) = c` requires `K_{i,c} = K_{j,pi(i)}` for
every `i` already assigned, which is `O(assigned)` per candidate and cuts at
once. In practice a positive answer costs about `N` nodes and a negative one
about `N^2`.

### Flags

| Flag | Meaning | Default |
|---|---|---|
| `--file FILE` | the matrix, as written by `get_chm`; **mandatory** | — |
| `--layout L` | `core` \| `full` \| `auto` | `core` |
| `--tol T` | two phases count as equal below `T` | `1e-9` |
| `--out FILE` | output path | `SH_<N>_<d>_<L>_<stamp>.m` |
| `--no-write` | decide only, write nothing | — |
| `--all` | try every base and report how many of the `N^2` work | — |
| `--quiet` | print the output file name only, nothing on failure | — |
| `--bare` | print `1` (symmetrisable) or `0`, and nothing else | — |
| `--help`, `-h` | usage and defaults | — |

### Exit codes

| Code | Meaning |
|---|---|
| 0 | a symmetric representative was found (and written) |
| 1 | the class contains none |
| 2 | bad usage, the file cannot be read, or it does not hold a CHM |

**A matrix that is not a CHM is refused, not answered.** Orthogonality is
used twice above — to know that column `b` is the only all-ones column of `K`,
which is what forces `pi(a) = b`, and to know that the witness is itself a CHM —
so on anything else both the verdict and the witness would be meaningless. The
usual cause is the wrong `--layout`, and the message says so.

**On the tolerance.** The two failure modes point in opposite directions: too
tight and distinct-but-close phases stop matching, giving a false NOT; too loose
and the codebook merges phases that differ, giving a false YES. The default
`1e-9` sits between the `1e-13` to which `get_chm` converges a matrix and the gap
between its distinct phases, which the run prints. When the gap approaches the
tolerance the run says so — at `N = 11` with 6051 distinct phases in the unit
interval a gap of `3e-8` is normal and still 30 times the tolerance, but it is
worth seeing. In the survey below every verdict was stable across `--tol` from
`1e-12` to `1e-3`.

```console
$ ./symmetrize --file YH_11_0_3081_20260905T074820.m
* YH_11_0_3081_20260905T074820.m   N = 11, nh = 6.247e-13, n1 = 1.570e-16
* codebook: 3081 distinct phases, smallest gap 7.005e-07 (tolerance 1e-09)
* base: dephase at (row 0, col 5)
* column permutation: 5 4 0 1 6 2 8 9 3 10 7
* || S - S^T ||_F = 8.277e-12,  nh = 6.242e-13,  n1 = 1.110e-16
* (11 search nodes)
* d = 0 and L = 3081 taken from the input file name (both are class invariants)
* SYMMETRISABLE   symmetric representative written to SH_11_0_3081_20260905T075508.m

$ ./symmetrize --file YH_11_0_6051_20260905T074716.m --bare
0
```

The column permutation is the whole witness: with `K` the matrix dephased at the
printed base, `S_ij = K_{i, pi(j)}` is the symmetric representative.

**Naming the result.** `d` and `#Lambda` are invariants of the equivalence
class, so when the input file is named in the `get_chm` convention and its `N`
agrees with the matrix actually read, the two numbers are carried over rather
than recomputed and only the leading letter changes to `S`. Otherwise `./defect`
and `./lambda` are run exactly as `get_chm` runs them. A hand-edited input name
is believed, so do not rename files by hand and then trust the output name.

**What a negative answer means.** "Not symmetrisable" is a statement about
*this* matrix. A matrix of positive defect sits in a continuous family whose other
members are inequivalent to it, and one of those may be symmetrisable when this
one is not — which is exactly what happens at `N = 6`, where the generic `d = 4`
member is not symmetrisable while the isolated `d = 0` Butson point is. For
defect 0 the class is isolated and the statement is about the class outright.

**Validation.** `check_me` confirms the necessary condition independently: a
symmetrisable matrix must be equivalent to its own transpose, and over 96
matrices at `N = 5 .. 11` the two verdicts agreed every time. A round trip also
holds — apply a random monomial transformation to a symmetric matrix, let
`symmetrize` recover a symmetric form, and `check_me` confirms it is the same
class, 18 of 18. The verdict is stable over `--tol` from `1e-12` to `1e-3`.

---

## Is every class symmetrisable?

**No.** The conjecture that for each `N` every equivalence class of CHM
contains a symmetric representative fails, and it fails already at small orders
and for isolated (defect 0) classes. Unconstrained searches:

| `N` | samples | symmetrisable | not |
|---|---|---|---|
| 5 | 12 | 12 | 0 |
| 6 | 12 | 1 | 11 |
| 7 | 12 | 12 | 0 |
| 8 | 12 | 2 | 10 |
| 9 | 12 | 8 | 4 |
| 10 | 12 | 3 | 9 |
| 11 | 24 | 3 | 21 |

At `N = 11` every sample had defect 0, so those 21 are isolated classes with no
symmetric representative at all. The strata behave uniformly: at `N = 11` the
`#Lambda = 3081` samples are symmetrisable and the `#Lambda = 6049`, `6051` and
`643` ones are not, with no stratum split between the two answers.

### An exact counterexample

One sample at `N = 9` is a `BH(9, 6)` — a Butson matrix over sixth roots of
unity, defect 0 — so it is described exactly by its integer LOG-form and no
floating point need enter the argument at all. Its LOG-form is

```
0 0 0 0 0 0 0 0 0
0 4 3 3 2 1 0 0 4
0 2 5 5 2 3 4 0 2
0 0 3 1 4 3 0 4 2      H_jk = exp(2 pi i * L_jk / 6)
0 2 4 2 4 0 4 2 0
0 0 1 3 0 3 4 2 4
0 4 3 5 0 5 2 2 2
0 4 0 2 2 4 2 4 0
0 2 1 5 4 1 2 4 4
```

and `symmetrize` exhausts every base:

```console
$ ./summary --file bh9.data --layout full --all
* this is a CHM (nf = 6.53855e-15, n1 = 4.57757e-16)
* this is a BH(9, 6)
* #L = 6
* ud = 0

$ ./symmetrize --file bh9.data --layout full --all --no-write
* codebook: 6 distinct phases, smallest gap 1.667e-01 (tolerance 1e-09)
* all 81 bases exhausted: 81 ruled out by the row/column
*   multisets alone, the rest by 0 backtracking nodes
* NOT SYMMETRISABLE  (no monomial transformation of this matrix is symmetric)
$ echo $?
1
```

The verdict is a theorem about this matrix rather than a numerical observation.
Its phases are the six sixth roots of unity, so the codebook gaps are `1/6`
while the decision tolerance is `1e-9` — eight orders of margin, and no
perturbation within reach of double precision can move a phase from one
codebook entry to another.

### The obstruction

If `S = D1 P1 H P2 D2` is symmetric then `S^T = S` lies in the class of `H^T` as
well as of `H`, so

```
symmetrisable  ==>  H is equivalent to H^T
```

and that necessary condition is what fails in every sample above: `check_me`
reports each non-symmetrisable matrix as inequivalent to its own transpose. The
converse is **not** a theorem, and the gap has a name.

Let the equivalence group act by `(M1, M2).H = M1 H M2`, and let

```
sigma(M1, M2) = (M2^T, M1^T) ,   so that   (g.H)^T = sigma(g).H^T .
```

If `H^T = h.H` then `H = sigma(h)h.H`, so `a := sigma(h)h` lies in `Aut(H)`. A
short computation shows that a symmetric representative exists exactly when `h`
can be chosen inside its coset `h.Aut(H)` with

```
sigma(h) h = 1 ,
```

i.e. when the element of the extended group that implements the transpose can be
taken to be an **involution**. In the language of the design-theory literature
`h` is a *duality* and such an `h` is a *polarity*: `H` is equivalent to `H^T`
iff a duality exists, and to a symmetric matrix iff a polarity exists. The
difference between the two is the whole story.

### Why the two coincide when `Aut(H)` is small

`Aut(H)` always contains the scalars `Z = {(lambda I, lambda^-1 I)}`. Suppose that
is all it contains. Then `a = sigma(h)h` is scalar, and writing `h = (M1, M2)`
the two components of `a` give `M2 = lambda^-1 M1^-T` and then
`lambda = lambda^-1`, so `lambda = +/-1`. The two cases are

```
lambda = +1 :  M2 =  M1^-T ,  and  S = M1 H  is SYMMETRIC
lambda = -1 :  M2 = -M1^-T ,  and  S = M1 H  is ANTISYMMETRIC
```

and the second is impossible for a complex Hadamard matrix: `M1 H` has the
entries of `H` permuted and rotated, so all of them are unimodular, and an
antisymmetric matrix has a zero diagonal. Hence

```
|Aut(H)/U(1)| = 1   and   H ~ H^T   ==>   H is symmetrisable.
```

This is the exact analogue of the Frobenius–Schur indicator, with the symplectic
alternative ruled out by the fact that no entry may vanish. It also explains the
coincidence observed above: a counterexample needs a nontrivial automorphism
group, and a matrix produced by a random search generally has none.
`check_me --count` prints `|Aut(H)/U(1)|`; across the 138 matrices surveyed here
it was 1 for every generic sample and larger only on the structured ones.

**Circulants are always symmetrisable.** If `H_jk = f(k-j)` then, with `P`
the reversal permutation `j -> -j`, one has `H^T = P H P` and `sigma(h)h = 1` for
`h = (P, P)`, so `P H` is symmetric outright. The same holds for the circulant
core families, which is why every circulant stratum found here
(`#Lambda = 5, 11, 63, 139, 161, 331` at `N = 11`) is symmetric or trivially
symmetrisable.

### The gap is real, and it is known — for real Hadamard matrices

Lin and Wallis (*Congr. Numer.* **85** (1991) 73–79) first exhibited a Hadamard
matrix equivalent to its transpose but not to a symmetric one. Valtonen's
classification (arXiv:2608.06920, 2026) settles the counts up to order 32:

| Order | Classes | Duality (`~ H^T`) | Polarity (symmetrisable) |
|---|---|---|---|
| 12 | 1 | 1 | 1 |
| 16 | 5 | 3 | 3 |
| 20 | 3 | 3 | 3 |
| 24 | 60 | 12 | 11 |
| 28 | 487 | 101 | 91 |
| 32 | 13 710 027 | ? | 3 231 |

So the first order at which the two columns differ is **24**, with
exactly one such class, and by order 28 there are ten. Those examples transfer
verbatim to the complex setting: if `H` is real with entries `+/-1` and
`D1 P1 H P2 D2` is also real, then all the phases are forced to be `+/-1` up to
one global scalar, so `Aut` and the dualities are the same group modulo scalars
and the polarity condition is unchanged. **A real Hadamard matrix is
symmetrisable over the complex monomial group exactly when it is symmetrisable
over signed permutations** — so the order-24 example is already a complex Hadamard matrix
that is equivalent to its transpose and to no symmetric matrix.

What is **not** settled is whether the gap occurs for a genuinely complex matrix,
with phases that are not roots of unity of order 2. Nothing here rules it out;
the place to look is matrices with large `|Aut(H)/U(1)|`, not the generic ones a
random search produces.

---

## Input files accepted

All five analysis programs read the same kind of file. The parser skips
`function ... / end`, ignores `%` and `#` comments, honours `...` continuations,
accepts `;` or a newline as the row separator and `,` or whitespace as the column
separator. A bare numeric table with no brackets — a plain `.data` file — is read
as well.

`--layout core` (the default) treats the array as the `(N-1)x(N-1)` core and
borders it with ones, so a 6×6 file describes a matrix of order 7.
`--layout full` treats it as the full `NxN` phase array. `--layout auto` tries
both and keeps whichever lands closer to a CHM.

If every entry is within `[-2pi, 2pi]` but some exceed 1 in modulus, the array is
taken to be in radians and divided by `2pi`, so a hand-made file in radians is
not silently misread.

---

## Tolerances

| Quantity | Flag | Default | Used by |
|---|---|---|---|
| convergence of a search | `--eps` | `1e-13` | `get_chm` |
| CHM test on `nf`, `n1` | `--eps` | `1e-7` | `summary` |
| rank of `R` | `--tol` | `1e-8` | `defect`, `summary` |
| uniqueness of `Lambda` | `--eps` / `--lambda-eps` | `1e-8` | `lambda`, `summary` |
| Butson test | `--butson-tol` | `1e-8` | `summary` |
| largest Butson exponent | `--qmax` | `10000` | `summary` |
| equality of two phases | `--tol` | `1e-9` | `check_me`, `symmetrize` |
| symmetry test in the name | not configurable | `1e-7` | `get_chm` |

**The Butson tolerance is `1e-8`.** A phase known to `1e-14` gives an error of
about `q * 1e-13` at the `q`-th power, which already exceeds `1e-10` around
`q = 1000`. A stricter threshold would therefore reject genuine Butson matrices
at the larger exponents `--qmax 10000` admits.

**The Butson test works on phases, not on powers.** Forming `M.^q` accumulates
round-off as `q` grows. Since
`|exp(2j*pi*q*theta) - 1| = 2|sin(pi*q*theta)|`, testing
`max_jk 2|sin(pi*q*theta_jk)| <= tol` costs one sine per entry and is exact to
within a single rounding at any `q`.

---

## Existence constraints

**A Hermitian CHM exists only for `N` even or `N` a perfect square.**

```
H = H'  and  H*H' = N*I  give  H^2 = N*I, so the eigenvalues are +/-sqrt(N)
and     tr H = sqrt(N) * (n_+ - n_-)
but     tr H = sum_j H_jj  is a sum of N real unimodular numbers, an integer
```

For `N` not a perfect square `sqrt(N)` is irrational, which forces `n_+ = n_-`
and hence `N` even. There is therefore nothing to look for at
`N = 3, 5, 7, 11, 13, 15`, where a search would otherwise run for ever.
`get_chm` refuses the request instead:

```console
$ ./get_chm --size 7 --symmetry h
./get_chm: no Hermitian CHM of order 7 exists.
  ...
  Use an even order or a perfect square.
$ echo $?
2
```

Symmetric CHM carry no such obstruction and are found for every order tried.

---

## Validation

The implementation is checked against results known in closed form. Each of
these can be reproduced with the shipped programs — `./defect`, `./lambda` and
`./summary` on a Fourier or Sylvester matrix written out as a phase file.

- `d(F_p) = 0` for `p` prime; `d(F_4) = 1`; `d(F_6) = 4`
- `d(F_{p^k}) = 1 + (k-1)*p^k - k*p^{k-1}`, verified for `N = 4, 8, 16, 32, 9, 27, 25, 49`
- `d(H_N) = (N-1)(N-2)/2` for the Sylvester real Hadamard matrices, verified for `N = 2, 4, 8, 16, 32`
- `#Lambda(F_p) = p` for `p` prime; `1` belongs to `Lambda`
- `F_N` is `BH(N, N)`; the Sylvester `H_4` is `BH(4, 2)`
- `SL(CHM) = 1` exactly; reshuffle and partial transpose are involutions
- `SL(H^G1) = SL(H^G2)`
- dephasing keeps the border of ones and preserves `H = H^T` and `H = H'`
- the core-phase file round-trips through disk without loss
- the polar factor of a unitary matrix is the matrix itself
- both methods produce genuine CHM, plain and with either symmetry

The two defect formulae are independent of each other and of this code, and both
are reproduced exactly, which exercises the whole path: building `R`, its SVD and
the rank threshold.

---

## Open questions

1. `rwcp` is exact and fast, but its landscape traps it past about `N = 10`. If
   larger orders matter for this method, a basin-hopping or parallel-tempering
   outer loop over the core phases is the next step.
2. For a Butson matrix the invariants are `q`-th roots of unity, so `Lambda`
   could be computed exactly in integer arithmetic modulo `q`, with no tolerance
   at all. That would make `#Lambda` exact rather than `eps`-dependent, and is
   worth doing if Butson cases come to dominate.
3. The defect is reported for the dephased orbit, `d(U) = (N-1)^2 - rank(R)`.
   The undephased defect is `D(U) = d(U) + 2N - 1`; it is not currently exposed,
   and would be a flag on `defect` if it were wanted.
