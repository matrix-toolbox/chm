#!/usr/bin/env python3
"""
build_matrices.py -- collect every matrix record the Catalog already holds
into one data file, matrices.js, which browse.html filters client-side.

2026-09-13  Claude Opus 5

Normally invoked through ./update_browse.sh, which also reports what changed.

Nothing is authored by hand: the records are read back out of index.html,
the appendices and the file names.  Re-run after editing any of them:

    python3 build_matrices.py [-o matrices.js] [--stats path/to/dL_statistic.md]

A record is one *observation*: a family probed at one parameter vector is
one record, a named matrix is one record.  So "d = 4 and #L < 100" answers
with the parameter vectors that realise it, not merely with family names.

Record keys are kept short because there are a few thousand of them:
  n  N            l  #Lambda       t  type flags   p  parameter vector
  d  defect       q  Butson q      nm name         f  .m file (repo path)
  s  source       u  page URL      c  comment      g  generic/underlined
  k  a related page (the catalogue entry, or the Appendix A section)
  a  appendix
"""
import argparse
import collections
import datetime
import html
import json
import re
from pathlib import Path

ROOT = Path(__file__).parent
BLOB = "https://github.com/matrix-toolbox/chm/blob/main/"

def detag(t):
    """Cell text without markup.  <br> separates two statements inside one
    cell, so it becomes "; " -- unless the text already punctuates itself."""
    t = re.sub(r"<br\s*/?>", "; ", t, flags=re.I)
    t = html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", t)))
    t = re.sub(r"\s*([;.:,])\s*;\s*", r"\1 ", t)      # "...;" + <br>  ->  "...; "
    return t.strip().strip(";").strip()


def rec(**kw):
    r = {"n": None, "d": None, "l": None, "q": None, "t": "", "nm": "",
         "p": None, "g": False, "s": "", "a": "", "f": "", "u": "", "k": "",
         "c": ""}
    r.update(kw)
    return r


def cells(row):
    return [detag(c) for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S)]


def row_links(row):
    hrefs = re.findall(r'href="([^"]+)"', row)
    f = next((h.replace(BLOB, "") for h in hrefs if h.startswith(BLOB)), "")
    u = next((h for h in hrefs if "catalogue/" in h), "")
    return f, u.replace("../", "")


# ------------------------------------------------------- the Catalog proper
LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def sigma_letters(tex):
    """'\\Sigma\\in\\{A,B,C,...,H\\}' -> ABCDEFGH ; '' if the name is plain"""
    m = re.search(r"\\Sigma\\in\\\{([^}]*)\\\}", tex)
    if not m:
        return []
    items = [x.strip() for x in m.group(1).split(",")]
    if "..." in items:                       # A,B,C,...,H  ->  A..H
        last = items[-1]
        return list(LETTERS[: LETTERS.index(last) + 1]) if last in LETTERS else []
    return [x for x in items if x in LETTERS]


def plain(tex):
    """$X_6^{(2)}$ -> X6 ;  $F_2\\otimes F_3$ -> F2xF3"""
    s = tex.strip("$ ")
    s = re.sub(r"\^\{?\((\d+)\)\}?", "", s)               # drop the superscript
    s = re.sub(r"\\Sigma.*", "", s)                        # C_{7\\Sigma : ...} -> C_{7
    s = s.replace(r"\otimes", "x")          # browse.html prints it back as \u2297
    s = re.sub(r"\\[a-zA-Z]+", "", s)
    return re.sub(r"[_{}()\\$\s:,]", "", s)


def family_dim(tex):
    m = re.search(r"\^\{?\((\d+)\)\}?", tex)
    return int(m.group(1)) if m else None


def defect_set(span):
    """'$\\{[\\underline{5}, 11], 13\\}$' -> ([5..11, 13], generic 5)"""
    body = re.sub(r"<[^>]+>", "", span)
    gen = [int(m.group(1)) for m in re.finditer(r"\\underline\{(\d+)\}", body)]
    body = re.sub(r"\\underline\{(\d+)\}", r"\1", body)
    body = re.sub(r"_\{[^}]*\}", "", body)
    ds = []
    for m in re.finditer(r"\[\s*(\d+)\s*,\s*(\d+)\s*\]", body):
        ds += list(range(int(m.group(1)), int(m.group(2)) + 1))
    body = re.sub(r"\[[^\]]*\]", "", body)
    ds += [int(x) for x in re.findall(r"\b\d+\b", body)]
    return sorted(set(ds)), (gen[0] if gen else None)


def mfile(*names):
    for n in names:
        if n and (ROOT / "CHM" / (n + ".m")).exists():
            return "CHM/%s.m" % n
    return ""


def catalog(path):
    src = re.sub(r"<!--.*?-->", "", path.read_text(encoding="utf-8"), flags=re.S)
    div = re.search(r'<div class="catalogue">(.*?)</div>', src, re.S)
    out = []
    for para in re.split(r"<hr\s*/?>", div.group(1)):
        for p in re.split(r"<p\b[^>]*>", para)[1:]:
            anchors = re.findall(r'href="catalogue/([^"]+)\.html">(.*?)</a>', p, re.S)
            if not anchors:
                continue
            span = re.search(r"<span[^>]*>(.*?)</span>", p, re.S)
            ds, gen = defect_set(span.group(1)) if span else ([], None)
            dl = re.search(r'href="(CHM_dL/index\.html#[^"]+)"', p)
            for i, (page, tex) in enumerate(anchors):
                N = int(page[:2]) if page[:2].isdigit() else None
                butson = page.endswith("B")
                fd = family_dim(tex)
                base_name = plain(tex)
                letters = sigma_letters(tex) or [""]
                t = ("B" if butson else "") + \
                    ("F" if fd else ("I" if fd == 0 else ""))
                note = []
                if fd:
                    note.append("family of dimension %d" % fd)
                # the defect set annotates the leading entry of the line
                mine = ds if i == 0 else []
                for L in letters:
                    nm = base_name + L
                    f = mfile(nm, "%s_%s" % (nm, fd) if fd else None)
                    common = dict(n=N, t=t, nm=nm, s="cat", a="", f=f,
                                  u="catalogue/%s.html" % page, c="; ".join(note),
                                  k=dl.group(1) if dl else "",
                                  q=int(page[2:4]) if butson and page[2:4].isdigit() else None)
                    if not mine:
                        out.append(rec(**common))
                    else:
                        for dv in mine:
                            out.append(rec(**common, d=dv, g=(dv == gen)))
    return out


# ---------------------------------------------------------------- Appendix A
def appendix_A(path):
    src = path.read_text(encoding="utf-8", errors="replace")
    # The BH(N, q) defect table is a grid of counts, not a list of matrices;
    # its rows would otherwise be read as entries named "2", "3", ...
    src = re.sub(r'<table id="BH_table">.*?</table>', "", src, flags=re.S)
    out = []
    heads = [(m.start(), detag(m.group(1)))
             for m in re.finditer(r"<h2[^>]*>(.*?)</h2>", src, re.S)]
    # An id belongs to the heading it sits in.  Keying by the heading's own
    # offset stops a section without an id from inheriting the previous one's
    # and pointing the reader at the wrong place.
    ids = {m.start(): m.group(1)
           for m in re.finditer(r'<h2[^>]*id="([^"]+)"', src)}

    for m in re.finditer(r"<pre[^>]*>(.*?)</pre>", src, re.S):
        title = next((t for p, t in reversed(heads) if p < m.start()), "")
        title = re.sub(r"[\u2191\u2193]\s*$", "", title).strip()
        fam = re.split(r"[(\s]", title, 1)[0] or title
        # Usually the dimension sits in the family name (F4, D8B, ...).  A few
        # headings carry it only in the argument -- "BH(10, 6)(p1)" -- so fall
        # back to the first number anywhere in the heading.
        mN = re.search(r"\d+", fam)
        if not mN:                      # "BH(10, 6)(p1)": keep the whole name
            mN = re.search(r"\d+", title)
            fam = re.sub(r"\s+", "", re.split(r"\)\(", title)[0] + ")") if "(" in title else title
        N = int(mN.group()) if mN else None
        if N is None:
            continue
        hp = next((p for p, _ in reversed(heads) if p < m.start()), None)
        anchor = ids.get(hp, "")
        url = "CHM_dL/index.html" + ("#" + anchor if anchor else "")
        # Most blocks end  ... d  #L, but a few carry a trailing q column, so
        # the last two numbers are not the invariants.  Read the header and use
        # the column positions it gives.
        col_d = col_l = None
        for raw in html.unescape(m.group(1)).splitlines():
            line = raw.strip()
            if not line or set(line) <= set("-| ") or "<" in line:
                continue
            f = line.split()
            if "d" in f and any(x.startswith("#") for x in f):
                col_d = f.index("d")
                col_l = next(i for i, x in enumerate(f) if x.startswith("#"))
                continue
            note, nm = "", re.search(r"\b(?:observed but )?not recorded\b", line)
            if nm:
                note, line = "not recorded", line[: nm.start()].strip()
                f = line.split()
            num = lambda x: x.rstrip("*").lstrip("-").isdigit()
            if (col_d is not None and not note and len(f) > col_l
                    and num(f[col_d]) and num(f[col_l])):
                ints = [f[col_d], f[col_l]]
                par = f[:col_d]
            else:
                ints = [x for x in f if num(x)]
                par = [] if note else f[: len(f) - 2]
            if len(ints) < 2:
                continue
            out.append(rec(n=N, d=int(ints[-2].rstrip("*")), l=int(ints[-1].rstrip("*")),
                           t="F", nm=fam, s="A", a="A", p=par or None, u=url, c=note,
                           g=("*" in line) or any(x.startswith("r") for x in par)))

    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", src, re.S):
        c = cells(row)
        if len(c) < 3 or not c[1].isdigit() or c[0] in ("", "H"):
            continue
        f, u = row_links(row)
        q = re.search(r"BH\(\s*(\d+),\s*(\d+)\s*\)", c[3] if len(c) > 3 else "")
        mn = re.search(r"\d+", c[0])
        if not mn:
            continue
        out.append(rec(n=int(mn.group()), d=int(c[1]),
                       l=int(c[2]) if c[2].isdigit() else None,
                       q=int(q.group(2)) if q else None, t="B" if q else "I",
                       nm=c[0], s="A", a="A", f=f, u="CHM_dL/index.html", k=u,
                       c=c[3] if len(c) > 3 else ""))
    return out


# ---------------------------------------------- Appendix C (symmetric etc.)
def appendix_C(path):
    src = path.read_text(encoding="utf-8", errors="replace")
    out, sym, N = [], "", None
    pat = r'id="(Symmetric|Hermitian)_CHM"|<h3[^>]*>N\s*=\s*(\d+)</h3>|<tr[^>]*>(.*?)</tr>'
    for m in re.finditer(pat, src, re.S):
        if m.group(1):
            sym = "S" if m.group(1) == "Symmetric" else "H"
        elif m.group(2):
            N = int(m.group(2))
        else:
            c = cells(m.group(3))
            if len(c) < 3 or not c[1].isdigit() or not c[2].isdigit():
                continue
            f, u = row_links(m.group(3))
            q = re.search(r"BH\(\s*\d+,\s*(\d+)\s*\)", " ".join(c[3:]))
            out.append(rec(n=N, d=int(c[1]), l=int(c[2]), t=sym + ("B" if q else ""),
                           q=int(q.group(1)) if q else None, nm=c[0], s="C", a="C", f=f,
                           u="CHM_SH/index.html#" + (sym.replace("S", "Symmetric")
                                                        .replace("H", "Hermitian") + "_CHM"),
                           k=u, c=" ".join(c[3:]).strip()))
    return out


# ------------------------------------------------------------- the .m files
# prefix -> type flags, per the schemes in the appendices
PREFIX = {"Y": "", "T": "", "xH": "", "BH": "B", "RH": "BR", "FH": "B",
          "SH": "S", "HH": "H", "TU": "U", "SR": "U", "SG": "U",
          "RD": "U", "GD": "U", "DS": "S", "YH": "H"}
KIND = {"TU": "2-unitary", "SR": "self R-dual", "SG": "self \u0393-dual",
        "RD": "R-dual", "GD": "\u0393-dual", "FH": "Fourier",
        "RH": "real Hadamard",
        "LH": "block-circulant, sequence L",
        "VH": "block-circulant, sequence V"}
DIR_APP = {"CHM_dL": "A", "CHM_SINKHORN": "B", "CHM_SH": "C",
           "CHM_kU": "D", "CHM_BC": "E", "CHM_BH_0": "A", "CHM": "",
           "CHM_GH": "F"}


def from_files():
    out = []
    # CHM_NUMERICAL/ is raw search output -- thousands of transient files that
    # come and go with every run, and are gitignored.  They are not catalog
    # entries; a matrix is kept by moving it into an appendix and renaming it.
    SKIP = {".git", "CHM_NUMERICAL"}
    for f in sorted(ROOT.rglob("*.m")):
        if SKIP & set(f.parts):
            continue
        m = re.fullmatch(r"([A-Za-z]+)_(\d+)_(\d+)_(\d+)(.*)", f.stem)
        if not m:
            continue
        pre, N, d, lam, suf = m.groups()
        d0 = f.relative_to(ROOT).parts[0]
        out.append(rec(n=int(N), d=int(d), l=int(lam), t=PREFIX.get(pre, ""),
                       nm=f.stem, s="file", a=DIR_APP.get(d0, ""),
                       f=str(f.relative_to(ROOT)),
                       u=(d0 + "/index.html") if (ROOT / d0 / "index.html").exists() else "",
                       c=KIND.get(pre, "")))

    # Appendix E keeps its matrices as .data.  The names there are generated, so
    # a strict match is safe; a .dat elsewhere carries a timestamp that would
    # be misread as an invariant (Y_11_0_7xx_20221115...).
    for f in sorted((ROOT / "CHM_BC").glob("*.data")):
        # Appendix E's names degrade as the invariants get too slow to compute:
        #   LH_N_d_L   both known      LH_N_d   #L unknown      LH_N   both unknown
        m = re.fullmatch(
            r"(LH|VH)_(\d+)(?:_(\d+)(?:_(\d+))?)?([A-Z])?(_[A-Za-z0-9_]+)?", f.stem)
        if not m:
            continue
        pre, N, d, lam, _letter, note = m.groups()
        c = KIND[pre]
        t = "C"
        if note:
            note = note[1:]
            c += "; " + note.replace("_", " ")
            if note.startswith("BH"):
                t += "B"
        out.append(rec(n=int(N), d=int(d) if d is not None else None,
                       l=int(lam) if lam is not None else None, t=t, nm=f.stem,
                       s="file", a="E", f=str(f.relative_to(ROOT)),
                       u="CHM_BC/index.html", c=c))

    # Appendix F keeps its matrices as .data too, named GH_N_d_L_group.  The
    # symmetry flag is not in the name, so it is read from GH_classes.tsv,
    # which gh_report.py writes alongside them.
    sym = {}
    tsv = ROOT / "CHM_GH" / "GH_classes.tsv"
    if tsv.exists():
        for line in tsv.read_text(encoding="utf-8").splitlines()[1:]:
            c_ = line.split("\t")
            if len(c_) >= 9 and c_[8]:
                sym[c_[8].rsplit(".", 1)[0]] = (c_[5] == "yes", c_[6] == "yes")
    # Appendix F stores generators, not matrices: one .gen row fixes H once the
    # group is known (see CHM_GH/GH_expand.m).
    for f in sorted((ROOT / "CHM_GH").glob("*.gen")):
        m = re.fullmatch(r"GH_(\d+)_(\d+)_(\d+)_(.+)", f.stem)
        if not m:
            continue
        N, d, lam, grp = m.groups()
        # GH_N_d_L_group__k when two classes share the whole label; the
        # trailing __k disambiguates the file, it is not part of the group.
        grp = re.sub(r"__\d+$", "", grp)
        sy, lat = sym.get(f.stem, (False, False))
        # L: the core is a Latin square of values, i.e. the generator is
        # injective -- no value repeats in any row or column of the core.
        t = "GI" + ("S" if sy else "") + ("L" if lat else "")
        out.append(rec(n=int(N), d=int(d), l=int(lam), t=t, nm=f.stem,
                       s="file", a="F", f=str(f.relative_to(ROOT)),
                       u="CHM_GH/index.html",
                       c="bordered group-developed over " + grp))
    return out


# ----------------------------------------------- numerically observed pairs
def from_stats(path):
    out, N, sym = [], None, ""
    for line in path.read_text(encoding="utf-8").splitlines():
        h = re.match(r"##\s*N\s*=\s*(\d+),\s*symmetry\s*=\s*(\w)", line)
        if h:
            N, sym = int(h.group(1)), h.group(2).upper()
            continue
        r = re.match(r"\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)", line)
        if r and N:
            out.append(rec(n=N, d=int(r.group(1)), l=int(r.group(2)),
                           t=sym if sym in ("S", "H") else "", s="num",
                           c="%s samples" % r.group(3)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", default="matrices.js")
    ap.add_argument("--stats", default="")
    a = ap.parse_args()

    recs = catalog(ROOT / "index.html")
    recs += appendix_A(ROOT / "CHM_dL" / "index.html")
    recs += appendix_C(ROOT / "CHM_SH" / "index.html")
    listed = {r["f"] for r in recs if r["f"]}
    files = from_files()
    bc = [r for r in files if r["a"] in ("E", "F")]  # E and F list them all
    orphans = [r for r in files if r["a"] not in ("E", "F")
               and r["f"] not in listed]
    for r in orphans:
        r["c"] = (r["c"] + "; " if r["c"] else "") + "not listed in the tables"
    recs += orphans + bc
    if a.stats and Path(a.stats).exists():
        recs += from_stats(Path(a.stats))

    by = {}
    for r in recs:
        by[r["s"]] = by.get(r["s"], 0) + 1
    print("  orphan .m files  :", len(orphans), "(named X_N_d_L, in no index)")
    print("  Appendix E .data :", sum(1 for r in bc if r["a"] == "E"))
    print("  Appendix F .data :", sum(1 for r in bc if r["a"] == "F"))
    print("  by source        :", ", ".join("%s=%d" % kv for kv in sorted(by.items())))
    print("  total            :", len(recs))

    print("  with parameters  :", sum(1 for r in recs if r["p"]))
    print("  d known          :", sum(1 for r in recs if r["d"] is not None))
    print("  #L known         :", sum(1 for r in recs if r["l"] is not None))
    print("  N range          :", min(r["n"] for r in recs if r["n"]),
          "..", max(r["n"] for r in recs if r["n"]))
    print("  distinct (N,d,#L):", len({(r["n"], r["d"], r["l"]) for r in recs}))

    for r in recs:                       # BH(N, 2) is a real Hadamard matrix
        if r["q"] == 2 and "R" not in r["t"]:
            r["t"] += "R"

    # d(H) = 0 implies H is isolated (Tadej-Zyczkowski).  The flag used to be
    # read off the catalogue superscript alone, so everything typed from a file
    # name -- the appendices and the orphans -- never got it.
    #
    # Only ever added, never removed: the converse does not hold.  A matrix can
    # be isolated with a non-zero defect -- the superscript (0) on C_6 and on
    # B_9 says exactly that, while their defects are 3 and 2.
    for r in recs:
        if r["d"] == 0 and "I" not in r["t"]:
            r["t"] += "I"

    # A name that resolves to a script in one place resolves to it everywhere.
    # The Catalog links A8 to CHM/A8.m; Appendix A lists the same matrix with no
    # link at all, so the row came out dead.  Only names with a single candidate
    # are filled in -- K6 points at both K6_2.m and K6_3.m and stays alone.
    cand = collections.defaultdict(set)
    for r in recs:
        if r["nm"] and r["f"]:
            cand[r["nm"]].add(r["f"])
    known = {n: v.pop() for n, v in cand.items() if len(v) == 1}
    filled = 0
    for r in recs:
        if r["nm"] and not r["f"] and r["nm"] in known:
            r["f"] = known[r["nm"]]
            filled += 1
    print("  links filled in  :", filled, "(name already resolved elsewhere)")
    print("  with .m file     :", sum(1 for r in recs if r["f"]))

    blank = rec()
    slim = [{k: v for k, v in r.items() if v != blank[k]} for r in recs]
    out = ROOT / a.o
    out.write_text('var matrices_built = "%s";\nvar matrices = %s;\n'
                   % (datetime.date.today().isoformat(),
                      json.dumps(slim, separators=(",", ":"))),
                   encoding="utf-8")
    print("  wrote %s (%d KB)" % (a.o, out.stat().st_size // 1024))


if __name__ == "__main__":
    main()
