#!/usr/bin/env python3
"""
link_sources.py -- link each catalogue heading to its Octave source

2026-09-13  Claude Opus 5

The Catalog lists its matrices but says nothing about the scripts that build
them.  This walks catalogue/*.html and, wherever a heading names a matrix for
which CHM/<name>.m exists, wraps that name in a link to the file on GitHub:

    <h1>$F_3^{(0)}$</h1>
    <h1><a href=".../blob/main/CHM/F3.m" ...>$F_3^{(0)}$</a></h1>

A heading may name several matrices ($C_{7A}$, $C_{7B}$, ...); each gets its
own link, and a Kronecker product links its factors ($F_2 \\otimes F_4$ leads to
F2.m and F4.m).  Headings whose matrix has no file are left alone, so the
script is safe to re-run after adding one:

    python3 link_sources.py [--dry]

It is idempotent: a heading that is already linked is skipped.
"""
import argparse
import os
import re

import build_matrices as bm

BLOB = "https://github.com/matrix-toolbox/chm/blob/main/"
MATH = re.compile(r"\$[^$]+\$")

# Names the rule cannot derive, each confirmed from the file's own header:
#   K9_2z.m  "Matrix K9_2 originally denoted as BC_9^{(2)}"
#   Q11X.m   "Symmetric isolated matrix Q11 ... d = 0 and #L = 63"
#   A15X.m   "Isolated CHM of order N = 15 found by A. Chan and A. Munemasa"
# and two given by W. Bruzda: V8A..V8D are solved in V8_ANALYTIC.m (V8X_0.m is
# marked obsolete there), and Y9 is the matrix stored as Y9,0,105.
ALIAS = {"K9": "CHM/K9_2z.m", "BC9": "CHM/K9_2z.m",
         "Q11A": "CHM/Q11X.m", "Q11B": "CHM/Q11X.m",
         "A15": "CHM/A15X.m",
         "V8A": "CHM/V8_ANALYTIC.m", "V8B": "CHM/V8_ANALYTIC.m",
         "V8C": "CHM/V8_ANALYTIC.m", "V8D": "CHM/V8_ANALYTIC.m",
         "Y9": "CHM_SINKHORN/Y_9_0_105.m"}


def source_of(seg):
    """The .m file a single $...$ heading fragment refers to, or ''."""
    name = bm.plain(seg)
    dim = bm.family_dim(seg)
    return bm.mfile(name, "%s_%s" % (name, dim) if dim else None) or \
        (ALIAS[name] if name in ALIAS and os.path.exists(ALIAS[name]) else "")


def anchor(f, math):
    return '<a href="%s%s" title="%s">%s</a>' % (BLOB, f, f, math)


def relink(head):
    """Wrap every named matrix in the heading; returns (html, [files]).

    A link may not sit inside $...$ -- it would break the typesetting -- so a
    Kronecker product is split into one math group per factor and the linked
    factors are joined by a bare $\\otimes$."""
    used = []

    def one(m):
        seg = m.group(0)
        f = source_of(seg)
        if f:
            used.append(f)
            return anchor(f, seg)
        if "\\otimes" not in seg:
            return seg
        out, hit = [], False
        for part in seg[1:-1].split("\\otimes"):
            math = "$" + part.strip() + "$"
            pf = source_of(math)
            if pf:
                used.append(pf)
                hit = True
                out.append(anchor(pf, math))
            else:
                out.append(math)
        # spaces matter: "$\\otimes$$H_2$" would read as a $$ display delimiter
        return " $\\otimes$ ".join(out) if hit else seg

    return MATH.sub(one, head), used


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry", action="store_true", help="report, change nothing")
    a = ap.parse_args()

    linked = plain = already = 0
    for name in sorted(os.listdir("catalogue")):
        if not name.endswith(".html") or name == "index.html":
            continue
        path = os.path.join("catalogue", name)
        src = open(path, encoding="utf-8").read()
        m = re.search(r"(<h1[^>]*>)(.*?)(</h1>)", src, re.S)
        if not m:
            continue
        if BLOB in m.group(2):
            already += 1
            continue
        head, used = relink(m.group(2))
        if not used:
            plain += 1
            continue
        linked += 1
        print("  %-16s -> %s" % (name, ", ".join(sorted(set(used)))))
        if not a.dry:
            open(path, "w", encoding="utf-8").write(
                src[:m.start()] + m.group(1) + head + m.group(3) + src[m.end():])

    print()
    print("  linked now       : %d" % linked)
    print("  already linked   : %d" % already)
    print("  no source file   : %d" % plain)
    if a.dry:
        print("  (dry run, nothing written)")


if __name__ == "__main__":
    main()
