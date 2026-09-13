# Keeping the Catalog up to date

A short working manual for <https://matrix-toolbox.github.io/chm>.

- [Where things are](#where-things-are)
- [Adding a matrix](#adding-a-matrix)
- [Updating the Browser](#updating-the-browser)
- [The numerical search](#the-numerical-search)
- [Git, the six commands you need](#git-the-six-commands-you-need)
- [Before you push](#before-you-push)

---

## Where things are

| | |
|---|---|
| `index.html` | the Catalog itself — the list of named matrices |
| `catalogue/NNxx.html` | one page per entry, 118 of them |
| `CHM/` | the Octave scripts those entries refer to |
| `chm_appendix.html` | the front page of the appendices |
| `CHM_dL/` | **Appendix A** — defect and #Λ |
| `CHM_SINKHORN/` | **Appendix B** — matrices found by Sinkhorn |
| `CHM_SH/` | **Appendix C** — symmetric and Hermitian |
| `CHM_kU/` | **Appendix D** — k-unitary, R-dual, Γ-dual |
| `CHM_BC/` | **Appendix E** — block-circulant |
| `CHM_BH_0/` | the Butson-Home defect tables |
| `CHM_tool/` | the C++ toolbox and its scripts — see `CHM_tool/CHM_SPECIFICATION.md` |
| `browse.html`, `matrices.js` | the Browser and its data — **generated, never edit by hand** |
| `literature.html`, `literature.js` | the bibliography; add a paper to the `.js` |

Each appendix directory has its own `index.html`. That page is the authority for
what the appendix contains; the Browser is built from it.

---

## Adding a matrix

**1. Name the file.** The convention across the whole Catalog is

```
<prefix>_<N>_<defect>_<#Lambda>.m          SH_9_0_681.m
```

with a trailing letter when two matrices share all three numbers
(`SH_10_0_349A`, `...B`) and a trailing word for a property worth recording
(`VH_11_0_161_persymmetric`). The prefixes in use:

```
SH symmetric    HH Hermitian    BH Butson       RH real Hadamard
LH, VH block-circulant (Appendix E, stored as .data)
TU 2-unitary    SR self R-dual  SG self Gamma-dual
Y, T, xH  found by Sinkhorn, no special symmetry
```

The **file name must equal the function name inside**, or Octave cannot call it.

**2. Get the invariants** — do not copy them from somewhere else.

The C++ tools read a **numeric phase array**, not an Octave function, so a
catalog `.m` file cannot be handed to them directly — it would fail with
*ragged array* or *no numeric data found*. Export the phases from Octave first:

```console
octave:1> addpath CHM_SH
octave:2> Y = SH_9_0_681;
octave:3> A = mod(angle(Y) / 2 / pi, 1);
octave:4> dlmwrite('/tmp/Y.data', A, 'delimiter', ' ', 'precision', 17)
```

```console
$ cd CHM_tool
$ ./summary --file /tmp/Y.data --layout full --all
* this is a CHM (nf = 1.32736e-14, n1 = 2.22045e-16)
* this is a symmetric matrix
* #L = 681
* ud = 0
```

`--layout full` because `A` above is the whole N×N array; the default `core`
expects the (N−1)×(N−1) core, which is what `get_chm` writes. Use `precision 17`
— `save -ascii` keeps only eight decimals and the residual then reads `1e-08`
instead of `1e-14`, which can shift the defect.

`./defect --file ... --bare` and `./lambda --file ... --bare --no-write` print
the two numbers alone, which is handy in a shell loop.

**3. Add a row** to the right appendix `index.html`. Copy a neighbouring row and
change it; every table has the same five columns:

```html
<tr><td><a href="https://github.com/matrix-toolbox/chm/blob/main/CHM_SH/SH_9_0_681.m">SH<sub>9,0,681</sub></a></td><td>0</td><td>681</td><td></td><td>comment</td></tr>
```

Rows are ordered by defect, then by #Λ.

**4. Rebuild the Browser** (next section) and commit.

---

## Updating the Browser

`browse.html` never changes. `matrices.js` is generated from `index.html`, the
appendix pages and the file names, so it goes stale the moment you edit any of
them. One command:

```console
$ ./update_browse.sh
  updated matrices.js -- reload browse.html
```

It reports what changed by name:

```
  new names    : SH9,0,18 SH9,0,309 SH9,0,681 SH9,16,3
```

and says `no change` when there is nothing to do, so it is safe to run any time.
To ask without writing anything:

```console
$ ./update_browse.sh -n
matrices.js is up to date.
```

That exits 1 when it is stale, which makes it a good last check before a commit.

**The one habit worth forming:** edit an index → run `./update_browse.sh`.

### Linking a heading to its script

`catalogue/NNxx.html` headings link to the Octave script that builds the matrix.
After adding a script to `CHM/`, run

```console
$ python3 link_sources.py           # --dry to see what it would do
```

It only ever adds a link when the file exists, and skips headings already
linked, so it is safe to re-run. Use `--force` after changing the alias table
inside it (the handful of names it cannot derive, such as `K9 -> K9_2z.m`).

---

## The numerical search

In `CHM_tool/`. Full details in `CHM_tool/CHM_SPECIFICATION.md`.

```console
$ ./search_chm.sh --count 10 --size 11 --symmetry s --lmin 1 --lmax 5700
```

Matrices land in `CHM_NUMERICAL/` beside `CHM_tool`, that is inside this
repository; one line per attempt is appended to `CHM_tool/CHM_NUMERICAL.tsv`.
**The matrices are disposable, the catalog is not** — it is the record of what
was searched, and the statistics are built from it alone. `CHM_NUMERICAL/` is
in `.gitignore`, so the raw output never lands in a commit; delete it whenever
it gets in the way. Keep a matrix by moving it into the appendix it belongs to
and giving it a proper name.

```console
$ ./show_histogram.sh          # CHM_NUMERICAL.tsv -> dL_statistic.md
$ ./update_SH_index.sh         # dL_statistic.md -> ../CHM_SH/index.html
```

`update_SH_index.sh` rewrites Appendix C in place, adding a greyed
*observed numerically* row for every (d, #Λ) pair the search found that the page
does not list yet. It never touches a pair that is already there, named or not,
and running it twice changes nothing. Then rebuild the Browser.

---

## Git, the six commands you need

```console
$ git status                 # what have I changed?
$ git diff                   # ... and exactly how?
$ git add -A                 # stage everything
$ git commit -m "short text" # record it, with a short subject line
$ git push                   # publish; the site updates a minute later
$ git pull                   # take in changes made elsewhere
```

Two more worth knowing:

```console
$ git checkout -- FILE       # throw away my changes to FILE
$ git log --oneline -5       # the last five commits
```

To stage only some of your changes, name them instead of `-A`:

```console
$ git add CHM_SH/index.html CHM_SH/SH_9_0_681.m matrices.js
```

Renaming a tracked file: use `git mv old new`, so the history follows it.

If `git push` is refused because the remote moved on, `git pull` first.

---

## Before you push

```console
$ ./update_browse.sh -n      # is matrices.js in sync?
$ git status                 # is anything unintended staged?
$ git diff --stat            # how big is this really?
```

Open `browse.html` and the page you edited in a browser. A page can be opened
straight from disk — the whole site works over `file://`, no server needed.

If a link in a table points at a file that does not exist, the Browser will show
the row but the link will 404 once pushed. A quick check that every `.m` or
`.data` referenced from an appendix is actually there:

```console
$ grep -o 'blob/main/[^"]*' CHM_SH/index.html | sed 's|blob/main/||' | \
      while read f; do [ -e "$f" ] || echo "missing: $f"; done
```

A blob link only resolves after the commit is **pushed** — a fresh file 404s on
GitHub until then, which is normal and not a broken link.

One trap worth knowing: this machine's locale writes a decimal **comma**, so
`awk`/`printf` in a shell pipeline can silently produce `0,5` where a tool
expects `0.5`. Octave's `dlmwrite` is not affected. Put `LC_ALL=C` in front of
any shell command that formats numbers — the scripts here already do.
