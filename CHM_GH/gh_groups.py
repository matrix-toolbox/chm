#!/usr/bin/env python3
"""
Complete library of all finite groups of order <= LIMIT, built from scratch.

Every group of order < 60 is solvable, hence has a normal subgroup of prime index.
So every G of order n arises as a cyclic extension

    G = H . <t>, H normal of index p, t^p = z in H, t h t^-1 = alpha(h)

for some group H of order n/p, some alpha in Aut(H) and some z in H with alpha(z) = z and alpha^p = inn_z.
Enumerating those pairs over all H of order n/p and all primes p | n is therefore *complete*.

Two reductions keep the count small:
  * alpha only matters up to Aut(H)-conjugacy  (beta: (alpha,z) -> (beta alpha beta^-1, beta(z)) is realised by h t^i -> beta(h) t^i);
  * isomorphic results are merged by an isomorphism-invariant fingerprint.

Since the fingerprint is a genuine invariant, distinct fingerprints always mean non-isomorphic groups,
so the count produced is a lower bound for the true number of groups of that order.
Agreement with the known counts (A000001), checked by self_test(), therefore certifies both completeness and correctness.
"""
import itertools, json, numpy as np, os, sys
from collections import Counter

# number of groups of order n, n = 1..34   (OEIS A000001)
A000001 = [0,1,1,1,2,1,2,1,5,2,2,1,5,1,2,1,14,1,5,1,5,2,2,1,15,2,2,5,4,1,4,1,51,1,2]


# ---------------------------------------------------------------- basics ----
def cyclic(n):
    i = np.arange(n)
    return (i[:, None] + i[None, :]) % n


def inverses(T):
    m = len(T)
    inv = np.empty(m, dtype=int)
    for i in range(m):
        inv[i] = int(np.nonzero(T[i] == 0)[0][0])
    return inv


def orders(T):
    m = len(T)
    o = np.empty(m, dtype=int)
    for i in range(m):
        k, x = 1, i
        while x != 0:
            x = T[x, i]
            k += 1
        o[i] = k
    return o


def is_group(T):
    m = len(T)
    if T.shape != (m, m):
        return False
    for r in T:                                   # latin square
        if len(set(r.tolist())) != m:
            return False
    for c in T.T:
        if len(set(c.tolist())) != m:
            return False
    if not (np.all(T[0] == np.arange(m)) and np.all(T[:, 0] == np.arange(m))):
        return False
    return bool(np.all(T[T, :][np.arange(m)] is not None) and
                np.all(T[T.ravel()][:, :].reshape(m, m, m) ==  # (ab)c
                       T[:, T].transpose(0, 1, 2)))            # a(bc)


def assoc_ok(T):
    m = len(T)
    left = T[T.reshape(-1), :].reshape(m, m, m)       # (a*b)*c
    right = T[:, T].reshape(m, m, m)                  # a*(b*c)
    return bool(np.array_equal(left, right))


# ------------------------------------------------------- generating sets ----
def closure(T, gens):
    seen = {0}
    frontier = [0]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = int(T[x, g])
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        frontier = nxt
    return seen


def min_gens(T):
    """A small generating set (greedy: always add the element that grows the
    generated subgroup the most)."""
    m = len(T)
    gens, cur = [], {0}
    while len(cur) < m:
        best, bestsz = None, len(cur)
        for g in range(m):
            if g in cur:
                continue
            sz = len(closure(T, gens + [g]))
            if sz > bestsz:
                best, bestsz = g, sz
        gens.append(best)
        cur = closure(T, gens)
    return gens


def word_table(T, gens):
    """BFS from the identity; return for each element a word (tuple of gen
    indices) spelling it."""
    m = len(T)
    word = [None] * m
    word[0] = ()
    frontier = [0]
    while frontier:
        nxt = []
        for x in frontier:
            for k, g in enumerate(gens):
                y = int(T[x, g])
                if word[y] is None:
                    word[y] = word[x] + (k,)
                    nxt.append(y)
        frontier = nxt
    return word


# ------------------------------------------------------------ Aut(H) -------
def automorphisms(T):
    """All automorphisms of T, each as a permutation array."""
    m = len(T)
    o = orders(T)
    gens = min_gens(T)
    words = word_table(T, gens)
    by_order = {}
    for x in range(m):
        by_order.setdefault(int(o[x]), []).append(x)

    out = []
    images = [None] * len(gens)

    def evaluate(imgs, k):
        """map determined by imgs[0..k-1]; returns array or None if inconsistent"""
        phi = np.full(m, -1, dtype=int)
        phi[0] = 0
        for x in range(m):
            w = words[x]
            if any(i >= k for i in w):
                continue
            y = 0
            for i in w:
                y = int(T[y, imgs[i]])
            if phi[x] != -1 and phi[x] != y:
                return None
            phi[x] = y
        return phi

    def rec(k, sub):
        if k == len(gens):
            phi = evaluate(images, k)
            if phi is None or len(set(phi.tolist())) != m:
                return
            # full homomorphism check
            if np.array_equal(phi[T], T[np.ix_(phi, phi)]):
                out.append(phi.copy())
            return
        for cand in by_order.get(int(o[gens[k]]), []):
            images[k] = cand
            newsub = closure(T, [images[i] for i in range(k + 1)])
            if len(newsub) != len(closure(T, gens[:k + 1])):
                continue
            phi = evaluate(images, k + 1)
            if phi is None:
                continue
            rec(k + 1, newsub)
    rec(0, {0})
    return out


def compose(a, b):
    return a[b]


def aut_class_reps(T, auts):
    """Representatives of the Aut(T)-conjugacy classes of Aut(T)."""
    key = {a.tobytes(): i for i, a in enumerate(auts)}
    inv = []
    for a in auts:
        ia = np.empty_like(a)
        ia[a] = np.arange(len(a))
        inv.append(ia)
    seen, reps = set(), []
    for i, a in enumerate(auts):
        if i in seen:
            continue
        reps.append(a)
        for j, b in enumerate(auts):
            c = compose(compose(b, a), inv[j])
            seen.add(key[c.tobytes()])
    return reps


# --------------------------------------------------------- the extension ---
def extend(T, p, auts=None):
    """All groups of order p*|T| having T as a normal subgroup of index p."""
    m = len(T)
    inv = inverses(T)
    if auts is None:
        auts = automorphisms(T)
    reps = aut_class_reps(T, auts)
    ident = np.arange(m)
    out = []
    for al in reps:
        # alpha^p
        ap = ident.copy()
        for _ in range(p):
            ap = compose(al, ap)
        fixed = [z for z in range(m) if al[z] == z]
        for z in fixed:
            innz = T[np.ix_([z] * m, ident)][0] if False else T[z][T[ident, inv[z]]]
            if not np.array_equal(ap, innz):
                continue
            # powers of alpha
            pw = [ident]
            for _ in range(p - 1):
                pw.append(compose(al, pw[-1]))
            M = p * m
            G = np.empty((M, M), dtype=int)
            for i in range(p):
                ai = pw[i]
                for j in range(p):
                    k = i + j
                    blk = T[np.ix_(ident, ai)]          # h * alpha^i(h')
                    if k >= p:
                        blk = T[blk, z]
                        k -= p
                    G[i * m:(i + 1) * m, j * m:(j + 1) * m] = blk + k * m
            out.append(G)
    return out


# ------------------------------------------------------------ invariants ---
def fingerprint(T):
    m = len(T)
    inv = inverses(T)
    o = orders(T)
    # conjugacy classes / centralizers
    conj = np.empty((m, m), dtype=int)
    for g in range(m):
        conj[g] = T[T[g], inv]                       # g * x * g^-1  (row g fixed)
    conj = np.empty((m, m), dtype=int)
    for g in range(m):
        conj[g] = T[T[np.full(m, g), np.arange(m)], inv[g]]
    cls, seen = [], set()
    for x in range(m):
        if x in seen:
            continue
        c = set(conj[:, x].tolist())
        seen |= c
        cls.append(c)
    centre = [x for x in range(m) if len(set(conj[:, x].tolist())) == 1]
    # derived subgroup
    comm = set()
    for a in range(m):
        for b in range(m):
            comm.add(int(T[inv[T[a, b]], T[b, a]]))
    D = closure(T, sorted(comm))
    # centraliser order profile per element
    prof = []
    for x in range(m):
        C = [g for g in range(m) if T[g, x] == T[x, g]]
        prof.append((int(o[x]), len(C), tuple(sorted(int(o[g]) for g in C))))
    # commutator-order profile per element
    cprof = []
    for a in range(m):
        cprof.append(tuple(sorted(int(o[T[inv[T[a, b]], T[b, a]]]) for b in range(m))))
    # power map
    pw = []
    for k in range(2, m + 1):
        s = set()
        for x in range(m):
            y, c = x, 1
            while c < k:
                y = int(T[y, x]); c += 1
            s.add(y)
        pw.append(len(s))
    return (m,
            tuple(sorted(o.tolist())),
            len(centre), len(D),
            tuple(sorted(int(o[g]) for g in centre)),
            tuple(sorted(int(o[g]) for g in D)),
            tuple(sorted(len(c) for c in cls)),
            tuple(sorted(prof)),
            tuple(sorted(cprof)),
            tuple(pw))



def sub_table(T, S):
    """Cayley table of the subgroup on the element set S (0 must be in S)."""
    S = sorted(S)
    idx = {g: i for i, g in enumerate(S)}
    k = len(S)
    return np.array([[idx[int(T[a, b])] for b in S] for a in S], dtype=int)


def index_p_subgroups(T, p):
    """All subgroups of index p that are kernels of a hom G -> Z_p.  For a
    p-group these are exactly the maximal subgroups."""
    m = len(T)
    gens = min_gens(T)
    words = word_table(T, gens)
    out = []
    for vals in itertools.product(range(p), repeat=len(gens)):
        if all(v == 0 for v in vals):
            continue
        phi = np.array([sum(vals[i] for i in words[x]) % p for x in range(m)])
        if not np.array_equal(phi[T], (phi[:, None] + phi[None, :]) % p):
            continue
        K = [x for x in range(m) if phi[x] == 0]
        if len(K) == m // p:
            out.append(K)
    # dedupe
    seen, uniq = set(), []
    for K in out:
        k = tuple(K)
        if k not in seen:
            seen.add(k)
            uniq.append(K)
    return uniq


def fingerprint2(T):
    """fingerprint(), refined by the isomorphism types of the maximal
    subgroups (one recursion level, which is enough here)."""
    base = fingerprint(T)
    m = len(T)
    primes = [q for q in range(2, m + 1) if m % q == 0 and
              all(q % r for r in range(2, q))]
    subs = []
    for q in primes:
        for K in index_p_subgroups(T, q):
            subs.append((q, fingerprint(sub_table(T, K))))
    return (base, tuple(sorted(map(repr, subs))))


def is_abelian(T):
    return bool(np.array_equal(T, T.T))


def abelian_type(T):
    """invariant factors, for naming"""
    m = len(T)
    o = orders(T)
    # decompose by brute force: greedily pick largest order not in span
    facs, cur = [], {0}
    while len(cur) < m:
        best = max((x for x in range(m) if x not in cur), key=lambda x: o[x])
        facs.append(int(o[best]))
        cur = closure(T, [])
        gens = []
        # recompute span of chosen generators
        chosen = []
        for f in facs:
            pass
        break
    # simpler: use the standard primary decomposition via counting
    from collections import Counter
    cnt = Counter(o.tolist())
    # determine type by matching against all partitions
    def types(n):
        # all abelian types of order n as sorted tuples of invariant factors
        def parts(k):
            if k == 1:
                yield ()
                return
            for d in range(2, k + 1):
                if k % d:
                    continue
                for rest in parts(k // d):
                    if not rest or d <= rest[0]:
                        yield (d,) + rest
        seen = set()
        for pp in parts(n):
            seen.add(tuple(sorted(pp)))
        return seen
    for t in sorted(types(m)):
        Tt = None
        for c in t:
            Tt = cyclic(c) if Tt is None else direct(Tt, cyclic(c))
        if Counter(orders(Tt).tolist()) == cnt:
            return t
    return None


def direct(A, B):
    a, b = len(A), len(B)
    M = a * b
    G = np.empty((M, M), dtype=int)
    for i in range(a):
        for j in range(b):
            for k in range(a):
                for l in range(b):
                    G[i * b + j, k * b + l] = A[i, k] * b + B[j, l]
    return G


# ------------------------------------------------------------ the library --
def build(limit=34, verbose=True, FP=None):
    FP = FP or fingerprint2
    lib = {1: [np.zeros((1, 1), dtype=int)]}
    autcache = {}
    for n in range(2, limit + 1):
        found, fps = [], set()
        primes = [p for p in range(2, n + 1) if n % p == 0 and
                  all(p % q for q in range(2, p))]
        for p in primes:
            h = n // p
            for hi, H in enumerate(lib[h]):
                key = (h, hi)
                if key not in autcache:
                    autcache[key] = automorphisms(H)
                for G in extend(H, p, autcache[key]):
                    f = FP(G)
                    if f in fps:
                        continue
                    fps.add(f)
                    found.append(G)
        lib[n] = found
        if verbose:
            exp = A000001[n] if n < len(A000001) else None
            flag = "ok" if exp == len(found) else f"MISMATCH (expected {exp})"
            print(f"  order {n:3d}: {len(found):3d} groups   {flag}", flush=True)
    return lib


LIB_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gh_lib.json")


def _save(lib, path=LIB_CACHE):
    """JSON, not pickle: the cache is data, and unpickling is code execution."""
    with open(path, "w") as fh:
        json.dump({str(n): [t.tolist() for t in ts] for n, ts in lib.items()}, fh)


def _load(path=LIB_CACHE):
    """Read the cache and check it really is a library of groups.

    Coercing arbitrary JSON to int and trusting the advertised order lets a
    malformed file through: [[0.5]] filed under order 34 truncates to [[0]] and
    is silently accepted as a group of order 34.  Everything is verified here
    instead, so a bad cache is rebuilt rather than used."""
    with open(path) as fh:
        raw = json.load(fh)
    lib = {}
    for key, tables in raw.items():
        n = int(key)
        out = []
        for t in tables:
            a = np.array(t)
            if a.dtype.kind not in "iu" and not np.all(a == np.floor(a)):
                raise ValueError(f"cache: non-integer table at order {n}")
            a = a.astype(int)
            if a.shape != (n, n):
                raise ValueError(f"cache: table {a.shape} filed under order {n}")
            if a.min() < 0 or a.max() >= n:
                raise ValueError(f"cache: entry out of range at order {n}")
            if not is_group(a):
                raise ValueError(f"cache: table at order {n} is not a group")
            out.append(a)
        lib[n] = out
    return lib


def load_or_build(limit=34, verbose=False):
    """The library, from the cache beside this file or freshly built into it.
    gh_sweep.py and gh_check.py both come through here, so the pipeline on the
    appendix page works from a clean checkout."""
    if os.path.exists(LIB_CACHE):
        try:
            lib = _load()
            if max(lib) >= limit:
                return lib
        except (ValueError, KeyError, json.JSONDecodeError):
            pass
    lib = build(limit, verbose=verbose)
    try:
        _save(lib)
    except OSError:
        pass
    return lib


if __name__ == "__main__":
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else 34
    print(f"building all groups of order <= {lim}")
    lib = build(lim)
    _save(lib)
    print(f"saved {sum(len(v) for v in lib.values())} groups to {LIB_CACHE}")
