#!/usr/bin/env python3
"""
Extended sweep for bordered group-developed complex Hadamard matrices.

For a finite group G of order m = N-1 and F : G -> T, the bordered matrix

    H = [ 1      1          ]
        [ 1   F(h g^-1)     ]_{g,h in G}

is complex Hadamard iff

    sum_g F(g) = -1
    A(v) := sum_k F(k) conj(F(kv)) = -1  (v != e).

There are m unknowns, m independent real conditions: critically determined for every G.
This script solves that system from many random starts,
for every group of order m and measures the resulting matrices with the catalogue tools.
"""
import numpy as np, sys, os, subprocess, json, hashlib
from scipy.optimize import least_squares

TOOL = os.environ.get("CHM_TOOL",
       os.path.join(os.path.dirname(os.path.dirname(
           os.path.abspath(__file__))), "CHM_tool"))


# ------------------------------------------------------------- the system --
class System:
    def __init__(self, T):
        self.T = T
        self.m = len(T)
        self.inv = np.array([int(np.nonzero(T[i] == 0)[0][0]) for i in range(self.m)])
        self.Tinv = T[:, self.inv]            # Tinv[m, v] = g_m * g_v^{-1}

    def A(self, F):
        return np.einsum('k,kv->v', F, np.conj(F[self.T]))

    def res(self, t):
        F = np.exp(1j * t)
        a = self.A(F)[1:] + 1.0
        s = F.sum() + 1.0
        return np.concatenate([[s.real, s.imag], a.real, a.imag])

    def jac(self, t):
        F = np.exp(1j * t)
        # dA(v)/dt_m = i[ F_m conj(F_{mv}) - F_{m v^-1} conj(F_m) ]
        J = 1j * (F[:, None] * np.conj(F[self.T]) -
                  F[self.Tinv] * np.conj(F)[:, None])       # (m unknown, v)
        J = J[:, 1:].T                                      # (v, m)
        g = 1j * F
        return np.vstack([[g.real], [g.imag], J.real, J.imag])

    def polish(self, t, rounds=80):
        """Gauss-Newton with a rank-revealing step: converges on the degenerate
        solution manifolds where Levenberg-Marquardt only crawls."""
        for _ in range(rounds):
            b = self.res(t)
            if np.abs(b).max() < 1e-15:
                break
            dt, *_ = np.linalg.lstsq(self.jac(t), -b, rcond=1e-10)
            t = t + dt
        return t

    def matrix(self, F):
        """the full N x N matrix, dephased"""
        m = self.m
        H = np.ones((m + 1, m + 1), dtype=complex)
        H[1:, 1:] = F[self.T[:, self.inv].T]     # H[i,j] = F(g_j g_i^-1)
        H = H / (H[:, :1] * H[:1, :] / H[0, 0])
        return H


def solve_all(T, restarts, seed, tol=1e-13):
    S = System(T)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(restarts):
        t0 = rng.uniform(0, 2 * np.pi, S.m)
        r = least_squares(S.res, t0, jac=S.jac, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        t = S.polish(r.x)
        if np.abs(S.res(t)).max() < tol:
            out.append(np.exp(1j * t))
    return S, out


# -------------------------------------------------------------- measuring --
def haagerup(H, tol=1e-8):
    """the Haagerup set {H_ij H_kl conj(H_il) conj(H_kj)}; its size is #Lambda.

    tol MUST match the catalogue: CHM_tool/lambda defaults to --eps 1e-8 and
    every #Lambda recorded elsewhere in the Catalog is taken at that value.  A
    coarser tolerance merges invariants that are genuinely distinct and reports
    #Lambda too low, which silently makes a matrix look like one the catalogue
    already holds."""
    L = np.einsum('ij,kl,il,kj->ijkl', H, H, np.conj(H), np.conj(H)).ravel()
    a = np.rint(L.real / tol).astype(np.int64)
    b = np.rint(L.imag / tol).astype(np.int64)
    return int(np.unique(a * np.int64(100000007) + b).size)


def orbit_key(S, F, tol=1e-6):
    """A cheap invariant of F under right translation F(g) -> F(gh) and the
    action of Aut(G).  Both permute the rows and columns of H, so solutions
    sharing this key are monomially equivalent and one may stand for the rest.

    Conjugation must NOT be folded in.  It is not part of the Catalog's
    equivalence relation, which is H -> D1 P1 H P2 D2, and a matrix here is
    routinely inequivalent to its own conjugate: check_me separates
    GH_11_0_139_D5, GH_11_0_139_Z2xZ5 and GH_12_0_463_Z11 from theirs, so
    identifying F with conj(F) would discard one of each such pair.

    The key is a rounded invariant, so two inequivalent solutions could still
    collide in one bucket.  That makes this a sampling heuristic, not a lossless
    reduction: the search is a sample of the solution set, never an enumeration
    of it."""
    n = int(round(1 / tol))
    P = F[:, None] * np.conj(F[S.T])                 # P[k, v] = F(k)conj(F(kv))
    q = np.rint(np.angle(P) / (2 * np.pi) * n).astype(np.int64) % n
    return tuple(sorted(tuple(sorted(q[:, v].tolist())) for v in range(S.m)))


def write_phases(H, path):
    P = (np.angle(H) / (2 * np.pi)) % 1.0
    with open(path, "w") as fh:
        for r in P:
            fh.write(" ".join(f"{q:.17g}" for q in r) + "\n")


def measure(H, path):
    """Defect and Butson type.

    The defect is taken from the singular values of the first-order system at a tolerance validated against 40-digit arithmetic,
    not from the catalogue binary: on a defective bordered matrix that system carries singular values near 1e-10
    and the binary's threshold sits above them, which inflates the defect.
    The error is one-sided, so d=0 verdicts agree either way.
    """
    import gh_report
    write_phases(H, path)
    d, _gap = gh_report.defect_svd(H, tol=1e-13)
    s = subprocess.run([f"{TOOL}/summary", "--file", path, "--layout", "full",
                        "--symmetry", "--butson"], capture_output=True,
                       text=True).stdout
    but = "no" if "not a BH" in s else ("yes" if "BH(" in s else "?")
    return int(d), but


# ------------------------------------------------------------------ main ---
def _task(arg):
    m, i, T, nm, R, seed, tmp, cap, cap2 = arg
    N = m + 1
    S, sols = solve_all(T, R, seed)
    thin = {}
    for F in sols:
        thin.setdefault(orbit_key(S, F), F)
    # Every candidate that survived the orbit thinning is kept.  #Lambda is NOT
    # an equivalence key -- inequivalent matrices over the same group routinely
    # share it (185 such triples in this appendix) -- so collapsing on it here
    # silently deleted distinguishable solutions.  Equivalence is decided later,
    # in gh_report.py, by an exhaustive monomial search.
    cands = []
    for F in list(thin.values())[:cap2]:
        H = S.matrix(F)
        cands.append((haagerup(H), H, F))
    path = os.path.join(tmp, f"gh_{m}_{i}.txt")
    recs = []
    for L, H, F in sorted(cands, key=lambda t: t[0]):
        d, but = measure(H, path)
        P = (np.angle(H) / (2 * np.pi)) % 1.0
        # Certificate: the smallest singular value of the Jacobian of the
        # Hadamard system.  Bounded away from zero means Newton converged
        # quadratically, so F is accurate to machine precision and the rank of
        # the first-order system -- hence the defect -- is unambiguous.  Near
        # zero means F has stalled at sqrt(eps) and the defect may be
        # under-stated; see the discussion of N=8 over Z_7.
        sv = np.linalg.svd(S.jac(np.angle(F)), compute_uv=False)
        cond = float(sv.min() / sv.max())
        rec = dict(d=d, L=L, butson=but, cond=cond,
                   degenerate=bool(cond < 1e-6),
                   sym=bool(np.abs(H - H.T).max() < 1e-9))
        if d == 0 and len([r for r in recs if r["d"] == 0]) < cap:
            rec["H"] = [[float(x) for x in r] for r in P]
        recs.append(rec)
    iso = sum(1 for r in recs if r["d"] == 0)
    deg = sum(1 for r in recs if r["degenerate"])
    return dict(N=N, m=m, gi=i, group=nm, restarts=R, hits=len(sols),
                orbits=len(thin), distinct=len(recs), isolated=iso,
                degenerate=deg, sols=recs)


def main():
    import argparse, multiprocessing as mp
    ap = argparse.ArgumentParser()
    ap.add_argument("--lib", default=None,
                    help="group library JSON; built on demand when omitted")
    ap.add_argument("--orders", default="4-34")
    ap.add_argument("--restarts", type=int, default=0)
    ap.add_argument("--seed", type=int, default=20260916)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--cap", type=int, default=60)
    ap.add_argument("--cap2", type=int, default=400)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tmp", default="/tmp")
    a = ap.parse_args()
    lo, hi = (int(x) for x in a.orders.split("-"))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    if a.lib:
        import gh_groups
        lib = gh_groups._load(a.lib)
    else:
        import gh_groups
        lib = gh_groups.load_or_build(hi)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import gh_names
    names = gh_names.names_for(lib)

    tasks = []
    for m in range(lo, hi + 1):
        R = a.restarts or (5000 if m <= 16 else 2500 if m <= 24 else
                          1000 if m <= 28 else 400)
        for i, T in enumerate(lib[m]):
            tasks.append((m, i, T, names[(m, i)], R, a.seed + 1000 * m + i,
                          a.tmp, a.cap, a.cap2))
    tasks.sort(key=lambda t: -t[0])          # big groups first: better packing
    rows = []
    with mp.Pool(a.jobs) as pool:
        for r in pool.imap_unordered(_task, tasks):
            rows.append(r)
            print(f"N={r['N']:3d} |G|={r['m']:3d} {r['group']:<24s} "
                  f"R={r['restarts']:5d} conv={r['hits']:5d} orb={r['orbits']:4d} "
                  f"distinct={r['distinct']:4d} isolated={r['isolated']:4d} "
                  f"degen={r['degenerate']:3d}"
                  + ("" if r["hits"] else "   NO SOLUTION"), flush=True)
            json.dump(sorted(rows, key=lambda x: (x["N"], x["gi"])),
                      open(a.out, "w"))
    json.dump(sorted(rows, key=lambda x: (x["N"], x["gi"])), open(a.out, "w"))


if __name__ == "__main__":
    main()
