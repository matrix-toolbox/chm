#!/usr/bin/env python3
"""
L309.py -- the symmetric CHM(9) with #Lambda = 309, built from thirteen numbers.

    python3 L309.py                # double precision, standard library only
    python3 L309.py --digits 60    # arbitrary precision (needs mpmath)

===============================================================================
WHAT THE MATRIX LOOKS LIKE
===============================================================================
Dephased at the right entry, the 8x8 core takes only 13 distinct values out of
64:

    [ 1  1  1  1  1  1  1  1  1 ]
    [ 1  a  b  c  d  e  f  g  h ]
    [ 1  b  a  h  g  f  e  d  c ]
    [ 1  c  h  i  j  c  h  j  k ]
    [ 1  d  g  j  l  g  d  m  j ]
    [ 1  e  f  c  g  a  b  d  h ]
    [ 1  f  e  h  d  b  a  g  c ]
    [ 1  g  d  j  m  d  g  l  j ]
    [ 1  h  c  k  j  h  c  j  i ]

===============================================================================
WHY THIRTEEN
===============================================================================
The core is symmetric, and its automorphism group (simultaneous row + column
permutations; found by brute force over all 8! = 40320 permutations) is a
KLEIN FOUR-GROUP of order 4, generated on the core indices 0..7 by

    (0 1)(2 7)(3 6)(4 5)      and      (0 4)(1 5)(3 6)

with third involution (0 5)(1 4)(2 7).  Note the group preserves the split
{0,1,4,5} | {2,3,6,7} and acts simply transitively on the first block -- a
regular Klein action, which is where the rigidity comes from.

The group <transpose, Aut> has order 8 and cuts the 64 core positions into 13
orbits (sizes 8,8,8,8,4,4,4,4,2,2,2,2,2), so 13 numbers fix the whole matrix
instead of the 36 a generic symmetric 8x8 core needs.

All nine #Lambda = 309 matrices found in a random search have exactly this
structure (13 values, |Aut| = 4), so it belongs to the stratum, not to one
lucky representative.

For comparison within N = 9:
    #Lambda = 681   -> |Aut| = 2, 20 free phases
    #Lambda = 309   -> |Aut| = 4, 13 free phases     <-- this file
and at N = 11 the circulant strata need only 5-6.

===============================================================================
THE EQUATIONS
===============================================================================
With conj(z) = 1/z (every entry is unimodular), orthogonality of rows p, q is

    sum_k H[p][k] / H[q][k] = 0        for every p < q

giving 72 real equations in the 13 unknowns -- overdetermined, consistent, and
with an isolated solution (the defect is 0).  Newton pins it to any precision.
===============================================================================
"""
import argparse
import cmath
import math
import sys

LETTERS = "abcdefghijklm"

# the 8x8 core pattern, one row per line
PATTERN = ["abcdefgh",
           "bahgfedc",
           "chijchjk",
           "dgjlgdmj",
           "efcgabdh",
           "fehdbagc",
           "gdjmdglj",
           "hckjhcji"]

# phases of the thirteen values, to 45 digits: a = exp(2 pi i * T[0]), etc.
# used only as the Newton starting point, so the script reproduces THIS matrix
T = ["0.037639104518053880001535411209507010008970883",   # a
     "0.736053342279209651063662604315843879985341365",   # b
     "0.484360779273302566813894120704199556148299864",   # c
     "0.363760368805647162497454900370433523236410251",   # d
     "0.260532071274409128675694569174041458743553369",   # e
     "0.325939265680887977672805647710128494062758870",   # f
     "0.768112169345912274349180611080320040265559082",   # g
     "0.778447081097861878782555735553843029974967806",   # h
     "0.119754239304489727542747771833392731829703995",   # i
     "0.218475702513797693114406289529606638541749709",   # j
     "0.668715596910520003246457383939440520564410638",   # k
     "0.830307262037130271323708334654579885702218507",   # l
     "0.604201565846324521220384891118851492399295266"]   # m

M = 8            # core size
N = 9            # matrix order
IDX = [[LETTERS.index(PATTERN[i][j]) for j in range(M)] for i in range(M)]
# the two generators of the automorphism group, on core indices
AUT = [[1, 0, 7, 6, 5, 4, 3, 2], [4, 5, 2, 6, 0, 1, 3, 7]]


def build_H(x, one):
    """The 9x9 matrix: a border of ones around the patterned core."""
    H = [[one for _ in range(N)] for _ in range(N)]
    for i in range(M):
        for j in range(M):
            H[i + 1][j + 1] = x[IDX[i][j]]
    return H


def residual(x, one):
    """sum_k H[p][k]/H[q][k] for every p < q -- zero iff H is a CHM."""
    H = build_H(x, one)
    out = []
    for p in range(N):
        for q in range(p + 1, N):
            s = one - one
            for k in range(N):
                s += H[p][k] / H[q][k]
            out.append(s)
    return out


def jacobian(x, one):
    """d(residual)/dx_s.  Each unknown occurs at several core positions, so the
    contributions accumulate."""
    H = build_H(x, one)
    rows = []
    for p in range(N):
        for q in range(p + 1, N):
            row = [one - one for _ in range(len(x))]
            for k in range(N):
                t = H[p][k] / H[q][k]
                if p >= 1 and k >= 1:
                    row[IDX[p - 1][k - 1]] += t / H[p][k]
                if q >= 1 and k >= 1:
                    row[IDX[q - 1][k - 1]] -= t / H[q][k]
            rows.append(row)
    return rows


def solve_normal(J, r, one):
    """Least squares through the normal equations: 72 equations, 13 unknowns."""
    n = len(J[0])
    A = [[one - one for _ in range(n)] for _ in range(n)]
    b = [one - one for _ in range(n)]
    for row, rv in zip(J, r):
        for i in range(n):
            if row[i] == 0:
                continue
            b[i] -= row[i].conjugate() * rv
            for j in range(n):
                A[i][j] += row[i].conjugate() * row[j]
    Mx = [A[i][:] + [b[i]] for i in range(n)]
    for c in range(n):
        piv = max(range(c, n), key=lambda rr: abs(Mx[rr][c]))
        if abs(Mx[piv][c]) == 0:
            raise ZeroDivisionError("singular normal matrix")
        Mx[c], Mx[piv] = Mx[piv], Mx[c]
        for rr in range(c + 1, n):
            fac = Mx[rr][c] / Mx[c][c]
            if fac != 0:
                for k in range(c, n + 1):
                    Mx[rr][k] -= fac * Mx[c][k]
    z = [None] * n
    for rr in reversed(range(n)):
        s = Mx[rr][n]
        for k in range(rr + 1, n):
            s -= Mx[rr][k] * z[k]
        z[rr] = s / Mx[rr][rr]
    return z


def newton(digits):
    if digits <= 15:
        one = complex(1, 0)
        x = [cmath.exp(2j * cmath.pi * float(t)) for t in T]
        tol = 1e-14
        rstr = lambda v, d: f"{float(v):.{d}e}"
        cstr = lambda z, d: f"{z.real:+.{d}f}{z.imag:+.{d}f}j"
    else:
        try:
            import mpmath as mp
        except ImportError:
            sys.exit("--digits > 15 needs mpmath:  pip install mpmath")
        mp.mp.dps = digits + 20
        one = mp.mpc(1, 0)
        x = [mp.expjpi(2 * mp.mpf(t)) for t in T]
        tol = mp.mpf(10) ** (-(digits + 5))
        rstr = lambda v, d: mp.nstr(v, d)
        cstr = lambda z, d: mp.nstr(z, d)
    for _ in range(200):
        r = residual(x, one)
        if max(abs(v) for v in r) < tol:
            break
        step = solve_normal(jacobian(x, one), r, one)
        x = [x[i] + step[i] for i in range(len(x))]
        # |x| = 1 is not part of the holomorphic system; project back so that
        # round-off cannot drift off the unit circle
        x = [v / abs(v) for v in x]
    return x, one, rstr, cstr


def lambda_count(H, tol=1e-8):
    n = len(H)
    vals = []
    for i in range(n):
        for k in range(i + 1, n):
            for j in range(n):
                for l in range(n):
                    z = H[i][j] * H[k][l] / (H[i][l] * H[k][j])
                    vals.append(math.atan2(complex(z).imag, complex(z).real)
                                / (2 * math.pi) % 1.0)
    vals.append(0.0)
    vals.sort()
    reps = [vals[0]]
    for v in vals[1:]:
        if v - reps[-1] >= tol:
            reps.append(v)
    if len(reps) > 1 and (1 - reps[-1]) + reps[0] < tol:
        reps.pop()
    return len(reps)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split("=====")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--digits", type=int, default=15)
    ap.add_argument("--show", type=int, default=12)
    a = ap.parse_args()

    print("the matrix, with placeholders  (thirteen distinct entries):\n")
    print("   [ " + "  ".join(["1"] * N) + " ]")
    for i in range(M):
        print("   [ 1  " + "  ".join(PATTERN[i][j] for j in range(M)) + " ]")
    print("\n   symmetric; automorphism group = Klein four-group generated by")
    print("   (0 1)(2 7)(3 6)(4 5) and (0 4)(1 5)(3 6) on the core indices.")

    x, one, rstr, cstr = newton(a.digits)
    print(f"\nsolving the system by Gauss-Newton to {a.digits} digits:\n")
    for s in range(len(x)):
        ph = math.atan2(complex(x[s]).imag, complex(x[s]).real) / (2 * math.pi) % 1.0
        print(f"   {LETTERS[s]} = exp(2*pi*i * {ph:.{min(a.show, 17)}f})"
              f"   = {cstr(x[s], a.show)}")

    H = build_H(x, one)
    resid = max(abs(sum(H[p][k] * H[q][k].conjugate() for k in range(N))
                    - (N if p == q else 0)) for p in range(N) for q in range(N))
    uni = max(abs(abs(H[i][j]) - 1) for i in range(N) for j in range(N))
    sym = max(abs(H[i][j] - H[j][i]) for i in range(N) for j in range(N))
    aut = max(abs(H[g[i] + 1][g[j] + 1] - H[i + 1][j + 1])
              for g in AUT for i in range(M) for j in range(M))
    print("\nchecks:")
    print(f"   max |(H H* - 9 I)_pq|                = {rstr(resid, 3)}")
    print(f"   max | |H_ij| - 1 |                   = {rstr(uni, 3)}")
    print(f"   max |H_ij - H_ji|         (symmetric)= {rstr(sym, 3)}")
    print(f"   max |C_(gi)(gj) - C_ij|   (Aut-inv.) = {rstr(aut, 3)}")
    print(f"   #Lambda                              = {lambda_count(H)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
