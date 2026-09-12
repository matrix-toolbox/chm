#!/usr/bin/env python3
"""
L139.py -- the symmetric CHM(11) with #Lambda = 139, built from six numbers.

    python3 L139.py                # double precision, no dependencies
    python3 L139.py --digits 60    # arbitrary precision (needs mpmath)
    python3 L139.py --minpoly      # also recover a and f from the degree-26
                                   # polynomial, to show how that fits in

===============================================================================
WHAT THE MATRIX LOOKS LIKE
===============================================================================
Only SIX distinct numbers a, b, c, d, e, f appear.  The dephased matrix is a
border of ones around a 10x10 circulant whose generating vector is palindromic,

    (a, b, c, d, e, f, e, d, c, b)

so the core is symmetric about BOTH diagonals (bisymmetric):

    [ 1 1 1 1 1 1 1 1 1 1 1 ]
    [ 1 a b c d e f e d c b ]
    [ 1 b a b c d e f e d c ]
    [ 1 c b a b c d e f e d ]
    [ 1 d c b a b c d e f e ]
    [ 1 e d c b a b c d e f ]
    [ 1 f e d c b a b c d e ]
    [ 1 e f e d c b a b c d ]
    [ 1 d e f e d c b a b c ]
    [ 1 c d e f e d c b a b ]
    [ 1 b c d e f e d c b a ]

===============================================================================
WHAT a, b, c, d, e, f ARE
===============================================================================
They are the six unimodular solutions of six explicit equations.  Writing
conj(z) = 1/z (legitimate because |z| = 1), row-orthogonality of H says

  (E0)  a + 2b + 2c + 2d + 2e + f                                       = -1
  (E1)  a/b + b/a + b/c + c/b + c/d + d/c + d/e + e/d + e/f + f/e       = -1
  (E2)  a/c + c/a + b/d + d/b + c/e + e/c + d/f + f/d + 2               = -1
  (E3)  a/d + d/a + b/c + c/b + b/e + e/b + c/f + f/c + d/e + e/d       = -1
  (E4)  a/e + e/a + b/d + d/b + b/f + f/b + c/e + e/c + 2               = -1
  (E5)  a/f + f/a + 2(b/e + e/b) + 2(c/d + d/c)                         = -1

E1..E5 are the autocorrelations A_1..A_5 of the generating vector; every one of
them is automatically REAL because the vector is palindromic.  Each grouped pair
z + 1/z is just 2*Re(z).

===============================================================================
ARE THE ENTRIES ROOTS OF ONE POLYNOMIAL?   No.
===============================================================================
This is the part that is easy to misread.  Eliminating variables from the
system above (Singular, 15 primes near 2^31, CRT + rational reconstruction)
gives a 0-dimensional ideal with vdim = 156 and these minimal-polynomial
degrees over Q:

    a, f          degree  78
    b, c, d, e    degree 140

The degree-26 polynomial reported earlier is the minimal polynomial of
    y = z + 1/z = 2*Re(z)
and only TWO of its 26 roots belong to this matrix: y = 2Re(a) and y = 2Re(f).
So the degree-26 polynomial does give a and f, via z = (y +- sqrt(y^2-4))/2
(run with --minpoly to see it) -- but it says nothing about b, c, d, e, and
even knowing all six minimal polynomials you would still have to work out which
root of each goes with which.  156/78 = 2, so a does not even determine the
rest by itself.

The honest "analytic form" is therefore the SYSTEM (E0..E5), not a root
extraction: six equations, six unknowns, an isolated solution.  This script
solves it by Newton's method, which converges quadratically to as many digits
as you ask for.
===============================================================================
"""
import argparse
import cmath
import sys

# ----------------------------------------------------------------------------
# Phases of the six entries, to 45 digits: a = exp(2 pi i * T[0]), etc.
# Used only as the starting point, so that Newton lands on THIS solution
# (the system has 16 unimodular solutions in 4 equivalence classes).
# ----------------------------------------------------------------------------
T = ["0.777013652127784570159804966544993144935821901",   # a
     "0.567459710463099252460896101474802918075716045",   # b
     "0.895898406314598549183621392322661135061602246",   # c
     "0.234518399045349256725006599363736330417841532",   # d
     "0.300642537831253596952207906255938855260599665",   # e
     "0.666365060256361202079656433069113772723603593"]   # f

NAMES = "abcdef"
# which of a..f sits at each position of the length-10 generating vector
GEN = [0, 1, 2, 3, 4, 5, 4, 3, 2, 1]
# positions of the generating vector driven by each unknown
SUP = [[i for i, g in enumerate(GEN) if g == j] for j in range(6)]


# ----------------------------------------------------------------------------
# the six equations, and their Jacobian
# ----------------------------------------------------------------------------
def equations(x, one):
    """F_t = 0 encodes (E1)..(E5) and (E0). x is [a,b,c,d,e,f]."""
    g = [x[i] for i in GEN]
    out = []
    for t in range(1, 6):
        s = one - one                                   # zero of the right type
        for m in range(10):
            s += g[m] / g[(m - t) % 10]
        out.append(s + one)
    s = one - one
    for m in range(10):
        s += g[m]
    out.append(s + one)
    return out


def jacobian(x, one):
    """d(equations)/dx, differentiated by hand (the system is holomorphic)."""
    g = [x[i] for i in GEN]
    J = [[one - one for _ in range(6)] for _ in range(6)]
    for t in range(1, 6):
        for j in range(6):
            s = one - one
            for m in SUP[j]:                            # numerator occurrences
                s += one / g[(m - t) % 10]
            for q in SUP[j]:                            # denominator occurrences
                s -= g[(q + t) % 10] / g[q] ** 2
            J[t - 1][j] = s
    for j in range(6):
        J[5][j] = one * len(SUP[j])
    return J


def solve_linear(A, b):
    """Gaussian elimination with partial pivoting -- 6x6, so no numpy needed."""
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) == 0:
            raise ZeroDivisionError("singular Jacobian")
        M[col], M[piv] = M[piv], M[col]
        for r in range(col + 1, n):
            fac = M[r][col] / M[col][col]
            if fac != 0:
                for k in range(col, n + 1):
                    M[r][k] -= fac * M[col][k]
    z = [None] * n
    for r in reversed(range(n)):
        s = M[r][n]
        for k in range(r + 1, n):
            s -= M[r][k] * z[k]
        z[r] = s / M[r][r]
    return z


def newton(digits):
    """Return [a,b,c,d,e,f] to `digits` correct digits."""
    if digits <= 15:
        one = complex(1, 0)
        x = [cmath.exp(2j * cmath.pi * float(t)) for t in T]
        tol = 1e-15
        absf = abs
        tostr = lambda z, d: f"{z.real:+.{d}f}{z.imag:+.{d}f}j"
        rstr = lambda v, d: f"{float(v):.{d}e}"
    else:
        try:
            import mpmath as mp
        except ImportError:
            sys.exit("--digits > 15 needs mpmath:  pip install mpmath")
        mp.mp.dps = digits + 20
        one = mp.mpc(1, 0)
        x = [mp.expjpi(2 * mp.mpf(t)) for t in T]
        tol = mp.mpf(10) ** (-(digits + 5))
        absf = abs
        tostr = lambda z, d: mp.nstr(z, d)
        rstr = lambda v, d: mp.nstr(v, d)
    for _ in range(200):
        F = equations(x, one)
        if max(absf(v) for v in F) < tol:
            break
        step = solve_linear(jacobian(x, one), [-v for v in F])
        x = [x[i] + step[i] for i in range(6)]
    return x, one, tostr, rstr


# ----------------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------------
def placeholder_matrix():
    rows = [["1"] * 11]
    for j in range(10):
        rows.append(["1"] + [NAMES[GEN[(k - j) % 10]] for k in range(10)])
    return rows


def build_H(x, one):
    g = [x[i] for i in GEN]
    H = [[one for _ in range(11)] for _ in range(11)]
    for j in range(10):
        for k in range(10):
            H[j + 1][k + 1] = g[(k - j) % 10]
    return H


def lambda_count(H, tol=1e-8):
    """#Lambda: distinct values of H_ij H_kl / (H_il H_kj)."""
    import math
    n = len(H)
    vals = []
    for i in range(n):
        for k in range(i + 1, n):
            for j in range(n):
                for l in range(n):
                    z = H[i][j] * H[k][l] / (H[i][l] * H[k][j])
                    p = math.atan2(complex(z).imag, complex(z).real) / (2 * math.pi)
                    vals.append(p % 1.0)
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
    ap = argparse.ArgumentParser(description=__doc__.split("=====")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--digits", type=int, default=15)
    ap.add_argument("--minpoly", action="store_true",
                    help="recover a and f from the degree-26 polynomial")
    ap.add_argument("--show", type=int, default=12, help="digits to print")
    a = ap.parse_args()

    print("the matrix, with placeholders  (only six distinct entries):\n")
    for r in placeholder_matrix():
        print("   [ " + "  ".join(f"{s:>1s}" for s in r) + " ]")

    x, one, tostr, rstr = newton(a.digits)
    print(f"\nsolving (E0)..(E5) by Newton to {a.digits} digits:\n")
    print(f"   {'':4s} {'modulus':>22s}   {'phase/2pi':>24s}")
    for j in range(6):
        z = x[j]
        import math
        ph = math.atan2(complex(z).imag, complex(z).real) / (2 * math.pi) % 1.0
        print(f"   {NAMES[j]} =  {rstr(abs(z), a.show):>22s}   "
              f"* exp(2*pi*i * {ph:.{min(a.show,17)}f})")
    print("\n   as complex numbers:")
    for j in range(6):
        print(f"   {NAMES[j]} = {tostr(x[j], a.show)}")

    H = build_H(x, one)
    resid = max(abs(sum(H[p][k] * H[q][k].conjugate() for k in range(11))
                    - (11 if p == q else 0)) for p in range(11) for q in range(11))
    sym = max(abs(H[i][j] - H[j][i]) for i in range(11) for j in range(11))
    uni = max(abs(abs(H[i][j]) - 1) for i in range(11) for j in range(11))
    # bisymmetry is a property of the 10x10 CORE: the border of ones is not
    # persymmetric (row 0 would have to match column 10)
    C = [[H[i + 1][j + 1] for j in range(10)] for i in range(10)]
    csym = max(abs(C[i][j] - C[j][i]) for i in range(10) for j in range(10))
    cper = max(abs(C[i][j] - C[9 - j][9 - i]) for i in range(10) for j in range(10))
    print("\nchecks:")
    print(f"   max |(H H* - 11 I)_pq|                  = {rstr(resid, 3)}")
    print(f"   max | |H_ij| - 1 |                      = {rstr(uni, 3)}")
    print(f"   max |H_ij - H_ji|            (symmetric)= {rstr(sym, 3)}")
    print(f"   core: max |C_ij - C_ji|      (symmetric)= {rstr(csym, 3)}")
    print(f"   core: max |C_ij - C_(9-j)(9-i)| (persym)= {rstr(cper, 3)}   -> bisymmetric")
    print(f"   #Lambda                                 = {lambda_count(H)}")

    if a.minpoly:
        print("\n--- where the degree-26 polynomial fits in ---")
        try:
            import mpmath as mp
        except ImportError:
            sys.exit("--minpoly needs mpmath")
        mp.mp.dps = 60
        P = MINPOLY26
        roots = mp.polyroots([mp.mpf(c) for c in P], maxsteps=500, extraprec=400)
        real = sorted(r.real for r in roots if abs(mp.im(r)) < mp.mpf(10) ** -25)
        print(f"   the degree-26 polynomial has {len(roots)} roots, "
              f"{len(real)} of them real")
        for nm, idx in (("a", 0), ("f", 5)):
            z = x[idx]
            y = mp.mpc(z) + 1 / mp.mpc(z)
            hit = min(real, key=lambda r: abs(r - y.real))
            print(f"   2*Re({nm}) = {mp.nstr(y.real, 25)}   nearest root "
                  f"{mp.nstr(hit, 25)}   diff {mp.nstr(abs(hit - y.real), 3)}")
        print("   so a and f come from this polynomial via z = (y +- sqrt(y^2-4))/2;")
        print("   b, c, d, e do NOT -- they have their own degree-140 polynomials.")


# the irreducible degree-26 minimal polynomial of y = z + 1/z for this family,
# ascending powers of y  (see EXACT_ALGEBRA.txt / y_minpoly.txt)
MINPOLY26 = [
    1296, 107136, -54158724,
    341782884, 108003677903, -1127688006500,
    -45391004059101, 59785524478804, 9528317590740053,
    83758633655905140, -1095203766827676215, -13724352040661927632,
    90049167440067093230, 803042391699799644576, -3707301230566177733210,
    -25735372609932345806576, 76199843375796469074226, 499546221040189249311568,
    -700429481907060356244746, -5677955651542023346380660, 528252036998369576970899,
    32422139783481434611103012, 28978081995971510284129511, -56570460150657248601111940,
    -90658708452903926132983527, -15283240909730206100358820, 17173646148407100390498325
]   # descending powers of y

if __name__ == "__main__":
    main()
