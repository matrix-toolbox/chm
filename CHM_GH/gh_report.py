#!/usr/bin/env python3
"""Turn GH_sweep.json into data files, equivalence certificates and LaTeX."""
import json, os, re, glob, csv, subprocess, itertools, sys
import numpy as np

ROOT = os.environ.get("CHM_ROOT",
       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TOOL = f"{ROOT}/CHM_tool"
HERE = os.path.dirname(os.path.abspath(__file__))


# ------------------------------------------------------- the catalogue -----
# Directories that must never contribute to the catalogue: the working area of
# this appendix itself (CHM_GH) and the gitignored research scratch space
# (CHM_WIP).  Both hold files whose names encode (N, d, #L) triples of
# *candidate* matrices, so counting them would let the sweep match its own
# output against itself.
SELF = ("CHM_GH", "CHM_WIP", "CHM_GD")


def catalogue():
    """(N, d, #Lambda) triples recorded anywhere in the published catalogue."""
    trip = set()
    src = open(f"{ROOT}/matrices.js").read()
    for r in json.loads(src[src.index("["):src.rindex("]") + 1]):
        # Appendix F is this sweep's own output; see SELF above.
        if r.get("a") == "F":
            continue
        if "d" in r and "l" in r:
            trip.add((r["n"], r["d"], r["l"]))
    with open(f"{TOOL}/CHM_NUMERICAL.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                trip.add((int(r["N"]), int(r["ud"]), int(r["#L"])))
            except (TypeError, ValueError):
                pass
    for f in glob.glob(f"{ROOT}/CHM_*/*.data") + glob.glob(f"{ROOT}/CHM_*/*.m"):
        if any(d in f for d in SELF):
            continue
        m = re.search(r"_(\d+)_(\d+)_(\d+)", os.path.basename(f))
        if m:
            trip.add(tuple(int(x) for x in m.groups()))
    return trip


def lambda_of(H, eps="1e-8", tmp="/tmp"):
    """#Lambda from CHM_tool/lambda -- the same binary, at the same eps, that
    every other appendix is measured with.  The Python routine in gh_sweep.py
    buckets the real and imaginary parts on a grid, which is a quantisation and
    not a distance test: it both splits values across a bin boundary and merges
    values further apart than eps.  It is fine for thinning a candidate list and
    must not be used for a recorded label."""
    D = H / (H[:, :1] * H[:1, :] / H[0, 0])
    core = (np.angle(D) / (2 * np.pi)) % 1.0
    path = os.path.join(tmp, "gh_lambda_%d.core" % os.getpid())
    np.savetxt(path, core[1:, 1:], fmt="%.17g")
    try:
        r = subprocess.run([f"{TOOL}/lambda", "--file", path, "--layout", "core",
                            "--eps", eps, "--no-write", "--bare"],
                           capture_output=True, text=True, timeout=600)
        return int(r.stdout.strip())
    finally:
        try: os.remove(path)
        except OSError: pass


# ------------------------------------------------------------ equivalence --
def core_file(P, path):
    with open(path, "w") as fh:
        for row in np.asarray(P)[1:, 1:]:
            fh.write(" ".join(f"{x:.17g}" for x in row) + "\n")


def equivalent(pa, pb, tol="1e-7"):
    """True / False / None, where None means the question was not answered.

    check_me exits 0 for equivalent, 1 for inequivalent and 2 for bad usage or
    an unreadable file.  Collapsing every non-zero status to False silently
    turns a crash, a missing binary or a bad argument into a confident
    "inequivalent", which is the one answer that manufactures new matrices."""
    try:
        r = subprocess.run([f"{TOOL}/check_me", pa, pb, "--quiet", "--tol", tol],
                           capture_output=True, text=True, timeout=1800)
    except (OSError, subprocess.TimeoutExpired) as e:
        print(f"    !! check_me did not run on {pa} {pb}: {e}")
        return None
    if r.returncode == 0:
        return True
    if r.returncode == 1:
        return False
    print(f"    !! check_me exited {r.returncode} on {pa} {pb}: "
          f"{(r.stderr or r.stdout).strip()[:120]}")
    return None


# ------------------------------------------------------------------ main ---
def main(src, outdir, tmp):
    rows = json.load(open(src))
    cat = catalogue()
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(tmp, exist_ok=True)

    # ---- collect every isolated solution that carries a matrix
    # Only solutions carrying the non-degeneracy certificate are admitted: at a
    # degenerate point Newton stalls at sqrt(eps), the vanishing singular values
    # of the first-order system are lifted to that level, and the defect can be
    # under-stated -- which would manufacture spurious isolated matrices.
    iso, dropped = [], 0
    for r in rows:
        for s in r["sols"]:
            if s["d"] != 0 or "H" not in s:
                continue
            if s.get("degenerate"):
                dropped += 1
                continue
            # Recompute #Lambda here rather than trusting the sweep's value:
            # gh_sweep.haagerup() buckets on a grid and is only a thinning key,
            # while every #Lambda recorded in the Catalog comes from
            # CHM_tool/lambda.  Taking the label from the matrix is what keeps
            # this appendix on the same footing as the others.
            H_ = np.exp(2j * np.pi * np.asarray(s["H"], dtype=float))
            iso.append(dict(N=r["N"], group=r["group"], gi=r["gi"],
                            L=lambda_of(H_), butson=s["butson"], sym=s["sym"],
                            P=s["H"]))
    print(f"  admitted {len(iso)} isolated solutions; "
          f"dropped {dropped} that failed the non-degeneracy certificate")
    # ---- de-duplicate by monomial equivalence, dimension by dimension
    classes = []
    undecided = []
    for N in sorted({x["N"] for x in iso}):
        here = [x for x in iso if x["N"] == N]
        reps = []
        for x in here:
            p = os.path.join(tmp, f"c_{N}_{len(reps)}_{len(classes)}.hp")
            core_file(x["P"], p)
            x["path"] = p
            hit = None
            for y in reps:
                # #Lambda is an exact equivalence invariant, so different
                # values rule equivalence out; both are measured here by
                # lambda_of at one tolerance, so the comparison is consistent.
                if y["L"] != x["L"]:
                    continue
                verdict = equivalent(x["path"], y["path"])
                if verdict is None:
                    undecided.append((x["path"], y["path"]))
                    continue
                if verdict:
                    hit = y
                    break
            if hit is None:
                x["groups"] = {x["group"]}
                reps.append(x)
            else:
                hit["groups"].add(x["group"])
        classes.extend(reps)
        print(f"  N={N:3d}: {len(here):4d} isolated solutions -> "
              f"{len(reps):3d} equivalence classes", flush=True)
    if undecided:
        # An unanswered comparison would otherwise become a silent "these are
        # different", i.e. a matrix announced as new on the strength of a
        # crashed subprocess.  Refuse to publish a classification built on that.
        raise RuntimeError(
            f"{len(undecided)} equivalence comparisons did not return a verdict; "
            f"first: {undecided[0]}.  Fix check_me before exporting.")

    # ---- independent defect check, with a separation certificate
    worst = float("inf")
    for c in classes:
        H = np.exp(2j * np.pi * np.asarray(c["P"]))
        c["d_svd"], c["gap"] = defect_svd(H)
        worst = min(worst, c["gap"])
    bad = [c for c in classes if c["d_svd"] != 0]
    print(f"  independent defect check: {len(classes)-len(bad)}/{len(classes)} "
          f"confirmed isolated, smallest separation {worst:.1e}")
    for c in bad:
        print(f"    !! N={c['N']} {c['group']} #L={c['L']}: "
              f"defect binary 0, SVD {c['d_svd']}")
    classes = [c for c in classes if c["d_svd"] == 0]

    # ---- write the new ones out
    # A triple match is a nomination, not a verdict: equal (N, d, #Lambda) does
    # not imply equivalence, and catalogue() collects labels taken under other
    # tolerance conventions (Appendix E states 1e-6, this appendix 1e-8).  A
    # class suppressed here is therefore only "already represented in the
    # Catalog by this invariant", which is what the page claims -- deciding it
    # would need the catalogued matrix itself and a check_me run against it.
    new = [c for c in classes if (c["N"], 0, c["L"]) not in cat]
    suppressed = [c for c in classes if (c["N"], 0, c["L"]) in cat]
    if suppressed:
        print(f"  {len(suppressed)} classes share a triple with a Catalog entry "
              f"and are not exported; equivalence to it was NOT tested")
    # A bordered group-developed matrix is fixed by the group and by F alone, so
    # only the generator is stored: row 0 of the core is F(h) for h in G, since
    # H[e,h] = F(h).  GH_groups.tsv carries the index tables needed to expand it.
    def sane(g):
        return (g.replace(":", "s").replace("(", "").replace(")", "")
                 .replace(",", "_").replace("#", "_"))
    used = {}
    for c in new:
        stem = f"GH_{c['N']}_0_{c['L']}_{sane(c['group'])}"
        # Two inequivalent classes can share N, d, #Lambda and the group; the
        # label is not a key.  Disambiguate instead of overwriting.
        used[stem] = used.get(stem, 0) + 1
        nm = stem if used[stem] == 1 else f"{stem}__{used[stem]}"
        with open(os.path.join(outdir, nm + ".gen"), "w") as fh:
            fh.write(" ".join(f"{x:.17g}" for x in c["P"][1][1:]) + "\n")
        c["file"] = nm + ".gen"
    write_group_tables({sane(c["group"]) for c in new}, outdir)
    json.dump([{k: v for k, v in c.items() if k not in ("P", "path", "groups")}
               | {"groups": sorted(c["groups"])} for c in classes],
              open(os.path.join(outdir, "GH_classes.json"), "w"), indent=1)
    print(f"\n{len(classes)} isolated classes, {len(new)} not in the catalogue")
    return rows, classes, new, cat




# ------------------------------------------------------------------ LaTeX --
def tex_escape(s):
    return (s.replace("x", "{\\times}").replace("Z", "\\Z_")
             .replace(":", "{\\rtimes}").replace("#", "\\#"))


def gname(s):
    """group name -> math mode"""
    if s.startswith("G("):
        return r"\mathrm{" + s.replace("#", "") + "}"
    out = s
    for a, b in (("Dih(", r"\mathrm{Dih}("), ("SL(2,3)", r"\mathrm{SL}(2,3)"),
                 ("D4oZ4", r"D_4\!\circ\!\Z_4")):
        out = out.replace(a, b)
    out = re.sub(r"Z(\d+)", r"\\Z_{\1}", out)
    out = re.sub(r"\bD(\d+)", r"D_{\1}", out)
    out = re.sub(r"\bQ(\d+)", r"Q_{\1}", out)
    out = re.sub(r"\bSD(\d+)", r"\\mathrm{SD}_{\1}", out)
    out = re.sub(r"\bM(\d+)", r"M_{\1}", out)
    out = re.sub(r"\bDic(\d+)", r"\\mathrm{Dic}_{\1}", out)
    out = re.sub(r"\bHe(\d+)", r"\\mathrm{He}_{\1}", out)
    out = out.replace("A4", "A_4").replace("S4", "S_4").replace("S3", "S_3")
    out = out.replace("x", r"{\times}").replace(":", r"{\rtimes}")
    return out


def emit(rows, classes, new, cat, path):
    from collections import defaultdict
    byN = defaultdict(list)
    for r in rows:
        byN[r["N"]].append(r)
    clsN, newN = defaultdict(list), defaultdict(list)
    for c in classes:
        clsN[c["N"]].append(c)
    for c in new:
        newN[c["N"]].append(c)

    # ---- Table 1: groups over which the system has no solution
    fails = [r for r in rows if r["hits"] == 0]
    t1 = [r"\begin{center}", r"\begin{tabular}{rll}", r"\toprule",
          r"$N$ & $|G|$ & groups admitting no solution\\", r"\midrule"]
    for N in sorted({f["N"] for f in fails}):
        gs = ", ".join("$" + gname(f["group"]) + "$"
                       for f in sorted(byN[N], key=lambda x: x["gi"])
                       if f["hits"] == 0)
        t1.append(f"{N} & {N-1} & {gs}\\\\")
    t1 += [r"\bottomrule", r"\end{tabular}", r"\end{center}"]

    # ---- Table 2: the census, one row per dimension
    t2 = [r"\begin{center}", r"\begin{tabular}{rrrrrr}", r"\toprule",
          r"$N$ & groups & solvable & isolated classes & new & $\Lam$ range\\",
          r"\midrule"]
    for N in sorted(byN):
        rs = byN[N]
        cs = clsN.get(N, [])
        rng = (f"${min(c['L'] for c in cs)}$--${max(c['L'] for c in cs)}$"
               if cs else "---")
        t2.append(f"{N} & {len(rs)} & {sum(1 for r in rs if r['hits'])} & "
                  f"{len(cs)} & {len(newN.get(N, []))} & {rng}\\\\")
    t2 += [r"\midrule",
           f"total & {len(rows)} & {sum(1 for r in rows if r['hits'])} & "
           f"{len(classes)} & {len(new)} & \\\\",
           r"\bottomrule", r"\end{tabular}", r"\end{center}"]

    # ---- Table 3: one concrete new matrix per dimension
    t3 = [r"\begin{center}", r"\begin{tabular}{rlrccl}", r"\toprule",
          r"$N$ & $G$ & $\Lam$ & Butson & symmetric & file\\", r"\midrule"]
    for N in sorted(newN):
        c = min(newN[N], key=lambda c: c["L"])
        t3.append(f"{N} & ${gname(c['group'])}$ & {c['L']} & "
                  f"{'yes' if c['butson'] == 'yes' else 'no'} & "
                  f"{'yes' if c['sym'] else 'no'} & "
                  f"{{\\tt {c['file'].replace('_', chr(92)+'_')[:-5]}}}\\\\")
    t3 += [r"\bottomrule", r"\end{tabular}", r"\end{center}"]

    open(path, "w").write("\n".join(t1) + "\n%%SPLIT%%\n" + "\n".join(t2)
                          + "\n%%SPLIT%%\n" + "\n".join(t3) + "\n")
    return len(fails)


def write_group_tables(groups, outdir):
    """I(i,j) = position of g_j g_i^-1, one row per group, 1-based for Octave."""
    import gh_groups as G, gh_names as NM
    lib = G.build(34, verbose=False)
    names = NM.names_for(lib)
    def sane(g):
        return (g.replace(":", "s").replace("(", "").replace(")", "")
                 .replace(",", "_").replace("#", "_"))
    idx = {sane(names[(n, i)]): (n, i) for n in lib for i in range(len(lib[n]))}
    with open(os.path.join(outdir, "GH_groups.tsv"), "w") as fh:
        fh.write("# group\torder\trow-major I(i,j) = index of g_j * g_i^-1, 1-based\n")
        for g in sorted(groups):
            n, i = idx[g]
            T = lib[n][i]
            inv = np.array([int(np.nonzero(T[k] == 0)[0][0]) for k in range(n)])
            I = (T[:, inv]).T + 1
            fh.write(f"{g}\t{n}\t" + " ".join(map(str, I.ravel())) + "\n")


# ----------------------------------------------- an independent defect ------
def defect_svd(H, tol=1e-13):
    r"""Defect computed directly from \eqref{eq:defdef}, together with a
    separation certificate.

    (H o R)H^dagger is Hermitian iff  sum_c H_ac conj(H_bc) (R_ac - R_bc) = 0
    for every a < b (the diagonal is automatically real).  That is N(N-1) real
    equations in the N^2 entries of R; the defect is the nullity minus 2N-1.
    Returns (defect, gap), where gap is the ratio of the smallest singular
    value counted as non-zero to the largest counted as zero: a large gap means
    the rank, hence the defect, is not an artefact of the threshold.
    """
    N = H.shape[0]
    rowsr, rowsi = [], []
    for a in range(N):
        for b in range(a + 1, N):
            w = H[a] * np.conj(H[b])              # coefficient of R_ac / R_bc
            M = np.zeros((N, N), dtype=complex)
            M[a] += w
            M[b] -= w
            rowsr.append(M.real.ravel())
            rowsi.append(M.imag.ravel())
    A = np.vstack(rowsr + rowsi)
    s = np.linalg.svd(A, compute_uv=False)
    s = s[s > 0]
    # Rank by the largest gap in the spectrum, not by a fixed threshold: on a
    # defective matrix the smallest non-zero singular value can sit near 1e-10,
    # only a few orders above machine noise, and any absolute cut near there
    # misclassifies it.  Only ratios below `floor` are considered, so that the
    # ordinary decay of an honest spectrum is not mistaken for the rank drop.
    cut = s.max() * tol
    rank = int((s > cut).sum())
    nz, zs = s[s > cut], s[s <= cut]
    gap = float(nz.min() / zs.max()) if len(zs) and len(nz) else float("inf")
    return N * N - rank - (2 * N - 1), gap


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
