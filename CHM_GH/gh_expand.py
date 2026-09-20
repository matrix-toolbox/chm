#!/usr/bin/env python3
"""
Rebuild a bordered group-developed CHM from its generator.

   python3 gh_expand.py GH_13_0_89_A4.gen # prints the N x N phase array

A .gen file holds one row: the phases of F : G -> U(1) in the order of the group elements.
With I(i,j) = position of g_j g_i^-1, tabulated in GH_groups.tsv,

    H[*, *] = H[*, h] = H[g, *] = 1
    H[g, h] = F(h g^-1) = F[I[g, h]],

so the matrix is just F read through I.
Phases are fractions of 2*pi, as everywhere in the catalog: H = exp(2*pi*i*A).
"""
import numpy as np, os, re, sys


def tables(path=None):
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "GH_groups.tsv")
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        g, n, body = line.rstrip("\n").split("\t")
        n = int(n)
        out[g] = np.fromstring(body, dtype=int, sep=" ").reshape(n, n) - 1
    return out


def expand(genfile, tabs=None):
    tabs = tabs if tabs is not None else tables()
    stem = os.path.basename(genfile).rsplit(".", 1)[0]
    # GH_N_d_L_group.gen, or GH_N_d_L_group__k.gen when two classes share the
    # whole (N, d, #Lambda, group) label -- the trailing __k is not part of the
    # group name.
    grp = re.sub(r"__\d+$", "", stem.split("_", 4)[4])
    a = np.loadtxt(genfile).ravel()
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{genfile}: generator has non-finite entries")
    I = tabs[grp]
    if len(a) != len(I):
        raise ValueError(f"{genfile}: {len(a)} phases but group {grp} has "
                         f"order {len(I)}")
    want = int(stem.split("_")[1])
    if want != len(I) + 1:
        raise ValueError(f"{genfile}: name says N={want} but group {grp} gives "
                         f"N={len(I) + 1}")
    F = np.exp(2j * np.pi * a)
    N = len(F) + 1
    H = np.ones((N, N), dtype=complex)
    H[1:, 1:] = F[I]
    return H


if __name__ == "__main__":
    H = expand(sys.argv[1])
    A = (np.angle(H) / (2 * np.pi)) % 1.0
    for row in A:
        print(" ".join(f"{x:.17g}" for x in row))
