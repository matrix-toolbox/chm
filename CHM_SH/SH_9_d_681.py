#!/usr/bin/env python3
"""
L681.py -- the symmetric CHM(9) with #Lambda = 681, built from twenty numbers.

    python3 L681.py                # double precision, standard library only
    python3 L681.py --digits 60    # arbitrary precision (needs mpmath)

===============================================================================
WHAT THE MATRIX LOOKS LIKE
===============================================================================
Dephased at the right entry, the 8x8 core takes only 20 distinct values out of
64, in this pattern:

    [ 1  1  1  1  1  1  1  1  1 ]
    [ 1  a  b  c  d  e  f  g  h ]
    [ 1  b  i  j  k  l  m  n  c ]
    [ 1  c  j  i  m  n  k  l  b ]
    [ 1  d  k  m  o  p  q  r  f ]
    [ 1  e  l  n  p  s  r  t  g ]
    [ 1  f  m  k  q  r  o  p  d ]
    [ 1  g  n  l  r  t  p  s  e ]
    [ 1  h  c  b  f  g  d  e  a ]

The pattern is symmetric, and it is also invariant under the involution

    sigma = (1 8)(2 3)(4 6)(5 7)          (rows and columns simultaneously)

acting on the eight core indices.  Checked by brute force over all 8! = 40320
permutations: the automorphism group of the core is EXACTLY of order 2, so
sigma is the whole story.  Burnside then gives the orbit count,

    (64 + 8 + 0 + 8) / 4 = 20

which is why 20 letters suffice instead of the 36 a generic symmetric 8x8 core
would need.  All eight #Lambda = 681 matrices found in a random search have the
same 20-value structure, so this is a property of the stratum, not an accident.

===============================================================================
HOW THIS DIFFERS FROM THE N = 11, #Lambda = 139 CASE
===============================================================================
For N = 11 the small-Lambda matrices were CIRCULANT, which cut the parameter
count to six and made exact elimination possible (an explicit degree-26
minimal polynomial).  Here that does not happen:

  * the #Lambda = 681 core is NOT circulant, and no dephasing makes it so;
  * the circulant families of order 9 give #Lambda in {19, 21, 27, 73, 76, 89,
    105, 201} -- never 681;
  * the reduction stops at 20 free phases.

Twenty unimodular unknowns means 40 variables in the elimination ideal, against
13 for the N = 11 case.  A Groebner basis at that size is not realistic, so this
script gives the exact SYSTEM and a solver, not a minimal polynomial.  The
matrix is still isolated (defect 0), so the solution below is rigid: Newton
pins it to as many digits as you ask for.

===============================================================================
THE EQUATIONS
===============================================================================
With H the 9x9 matrix above and conj(z) = 1/z (all entries are unimodular),
orthogonality of rows p and q reads

    sum_k H[p][k] / H[q][k] = 0     for every p < q

which, after substituting the pattern, is a system in the twenty unknowns
a, b, ..., t.  It is overdetermined (72 real equations, 20 unknowns) and
consistent, which is exactly what an isolated solution of a symmetric ansatz
looks like.
===============================================================================
"""
import argparse
import cmath
import sys

LETTERS = "abcdefghijklmnopqrst"

# the 8x8 core pattern, one row per line
PATTERN = ["abcdefgh",
           "bijklmnc",
           "cjimnklb",
           "dkmopqrf",
           "elnpsrtg",
           "fmkqropd",
           "gnlrtpse",
           "hcbfgdea"]

# phases of the twenty values, to 45 digits: a = exp(2 pi i * T[0]), etc.
# used only as the Newton starting point, so the script reproduces THIS matrix
T = ["0.913346652444725342467192362433566395976443023",   # a
     "0.670772659169588712723295149935416963769776542",   # b
     "0.367291105158288733798851206974233421541633555",   # c
     "0.619552200613313476078480961940689966169478717",   # d
     "0.221606545678372835478341299757154035375023812",   # e
     "0.945387166744969340179978712400246846147708166",   # f
     "0.272350270386420640222938528851464615045201791",   # g
     "0.548370928426654463970281503365356855783740411",   # h
     "0.873589548462712068770013226954115578772301551",   # i
     "0.624678858896765049249067533367219708431122270",   # j
     "0.209041559429553615367162754629606361882251894",   # k
     "0.118519123752116893746342923051676886050937179",   # l
     "0.335432263983712802491530706928658275563228168",   # m
     "0.697958102729687504601980605652933168182087127",   # n
     "0.475293336491417024951215132062250358529539001",   # o
     "0.534722763630793172974336234706555681782997072",   # p
     "0.040778216017792328015367516080496159006962506",   # q
     "0.756839967337246004150506238147154476780008968",   # r
     "0.812450743591537785644721113449589035506680820",   # s
     "0.428460416486727919980902092915485468040238077"]   # t

M = 8            # core size
N = 9            # matrix order
IDX = [[LETTERS.index(PATTERN[i][j]) for j in range(M)] for i in range(M)]


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
    """d(residual)/dx_s, by hand.  Entry (i,j) of the core belongs to unknown
    IDX[i][j], so each unknown occurs in several places of every row."""
    H = build_H(x, one)
    rows = []
    for p in range(N):
        for q in range(p + 1, N):
            row = [one - one for _ in range(len(x))]
            for k in range(N):
                t = H[p][k] / H[q][k]
                if p >= 1 and k >= 1:
                    row[IDX[p - 1][k - 1]] += t / H[p][k]      # numerator
                if q >= 1 and k >= 1:
                    row[IDX[q - 1][k - 1]] -= t / H[q][k]      # denominator
            rows.append(row)
    return rows


def solve_normal(J, r, one):
    """Least squares via the normal equations: 72 equations, 20 unknowns."""
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
        # project back onto the unit circle; the constraint |x|=1 is not part
        # of the holomorphic system and round-off would slowly drift off it
        x = [v / abs(v) for v in x]
    return x, one, rstr, cstr


def lambda_count(H, tol=1e-8):
    import math
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

    print("the matrix, with placeholders  (twenty distinct entries):\n")
    print("   [ " + "  ".join(["1"] * N) + " ]")
    for i in range(M):
        print("   [ 1  " + "  ".join(PATTERN[i][j] for j in range(M)) + " ]")
    print("\n   symmetric, and invariant under sigma = (1 8)(2 3)(4 6)(5 7)")
    print("   acting on rows and columns of the core simultaneously.")

    x, one, rstr, cstr = newton(a.digits)
    print(f"\nsolving the system by Gauss-Newton to {a.digits} digits:\n")
    import math
    for s in range(len(x)):
        ph = math.atan2(complex(x[s]).imag, complex(x[s]).real) / (2 * math.pi) % 1.0
        print(f"   {LETTERS[s]} = exp(2*pi*i * {ph:.{min(a.show, 17)}f})"
              f"   = {cstr(x[s], a.show)}")

    H = build_H(x, one)
    resid = max(abs(sum(H[p][k] * H[q][k].conjugate() for k in range(N))
                    - (N if p == q else 0)) for p in range(N) for q in range(N))
    uni = max(abs(abs(H[i][j]) - 1) for i in range(N) for j in range(N))
    sym = max(abs(H[i][j] - H[j][i]) for i in range(N) for j in range(N))
    sig = [7, 2, 1, 5, 6, 3, 4, 0]
    aut = max(abs(H[sig[i] + 1][sig[j] + 1] - H[i + 1][j + 1])
              for i in range(M) for j in range(M))
    print("\nchecks:")
    print(f"   max |(H H* - 9 I)_pq|                = {rstr(resid, 3)}")
    print(f"   max | |H_ij| - 1 |                   = {rstr(uni, 3)}")
    print(f"   max |H_ij - H_ji|         (symmetric)= {rstr(sym, 3)}")
    print(f"   max |C_(si)(sj) - C_ij|  (sigma-inv.)= {rstr(aut, 3)}")
    print(f"   #Lambda                              = {lambda_count(H)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
