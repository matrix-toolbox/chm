#!/usr/bin/env python3
"""Readable names for the groups produced by gh_groups.build()."""
import numpy as np, itertools
from gh_groups import (cyclic, direct, orders, inverses, closure, is_abelian,
                       fingerprint2, sub_table, min_gens)


def _abelian_type(T):
    m = len(T)
    if m == 1:
        return ()
    from collections import Counter
    cnt = Counter(orders(T).tolist())

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
    for pp in sorted({tuple(sorted(x)) for x in parts(m)}):
        Tt = None
        for c in pp:
            Tt = cyclic(c) if Tt is None else direct(Tt, cyclic(c))
        if Counter(orders(Tt).tolist()) == cnt:
            return pp
    return None


def _fmt_abelian(t):
    return "x".join(f"Z{c}" for c in t) if t else "1"


def _dihedral(n):
    """D_n of order 2n: rotations r^i and reflections s r^i."""
    m = 2 * n
    T = np.empty((m, m), dtype=int)
    for a in range(m):
        for b in range(m):
            ia, sa = a % n, a // n
            ib, sb = b % n, b // n
            i = (ia + ib) % n if sa == 0 else (ia - ib) % n
            T[a, b] = i + ((sa + sb) % 2) * n
    return T


def _dicyclic(n):
    """Dic_n = Q_{4n}: <a,b | a^{2n}=1, b^2=a^n, b a b^-1 = a^-1>, order 4n."""
    m = 4 * n
    T = np.empty((m, m), dtype=int)
    for a in range(m):
        for b in range(m):
            ia, sa = a % (2 * n), a // (2 * n)
            ib, sb = b % (2 * n), b // (2 * n)
            i = (ia + ib) % (2 * n) if sa == 0 else (ia - ib) % (2 * n)
            s = sa + sb
            if s == 2:
                i = (i + n) % (2 * n)
            T[a, b] = i + (s % 2) * (2 * n)
    return T


def _semidihedral(k):
    """SD_{2^k}: <a,b | a^{2^{k-1}}=b^2=1, bab^-1 = a^{2^{k-2}-1}>"""
    n = 2 ** (k - 1)
    t = 2 ** (k - 2) - 1
    m = 2 * n
    T = np.empty((m, m), dtype=int)
    for a in range(m):
        for b in range(m):
            ia, sa = a % n, a // n
            ib, sb = b % n, b // n
            i = (ia + ib) % n if sa == 0 else (ia + t * ib) % n
            T[a, b] = i + ((sa + sb) % 2) * n
    return T


def _modular(k):
    """M_{2^k}: <a,b | a^{2^{k-1}}=b^2=1, bab^-1 = a^{2^{k-2}+1}>"""
    n = 2 ** (k - 1)
    t = 2 ** (k - 2) + 1
    m = 2 * n
    T = np.empty((m, m), dtype=int)
    for a in range(m):
        for b in range(m):
            ia, sa = a % n, a // n
            ib, sb = b % n, b // n
            i = (ia + ib) % n if sa == 0 else (ia + t * ib) % n
            T[a, b] = i + ((sa + sb) % 2) * n
    return T


def _semidirect_cyclic(a, b, t):
    """Z_a : Z_b with generator of Z_b acting as x -> x^t."""
    m = a * b
    T = np.empty((m, m), dtype=int)
    for x in range(m):
        for y in range(m):
            ix, jx = x % a, x // a
            iy, jy = y % a, y // a
            T[x, y] = ((ix + pow(t, jx, a) * iy) % a) + ((jx + jy) % b) * a
    return T


def _A4():
    perms = [p for p in itertools.permutations(range(4))
             if sum(1 for i in range(4) for j in range(i) if p[i] < p[j]) % 2 == 0]
    idx = {p: i for i, p in enumerate(perms)}
    perms.sort(key=lambda p: 0 if p == (0, 1, 2, 3) else 1)
    idx = {p: i for i, p in enumerate(perms)}
    m = len(perms)
    return np.array([[idx[tuple(perms[a][perms[b][i]] for i in range(4))]
                      for b in range(m)] for a in range(m)], dtype=int)


def _perm_group(perms):
    perms = sorted(perms, key=lambda p: (p != tuple(range(len(p))), p))
    idx = {p: i for i, p in enumerate(perms)}
    m = len(perms)
    return np.array([[idx[tuple(perms[a][perms[b][i]] for i in range(len(perms[0])))]
                      for b in range(m)] for a in range(m)], dtype=int)


def _Sn(n):
    return _perm_group(list(itertools.permutations(range(n))))


def _SL23():
    """SL(2,3) as 2x2 matrices over F_3 of determinant 1, order 24."""
    els = [(a, b, c, d) for a in range(3) for b in range(3)
           for c in range(3) for d in range(3) if (a * d - b * c) % 3 == 1]
    idx = {e: i for i, e in enumerate(sorted(els, key=lambda e: e != (1, 0, 0, 1)))}
    els = sorted(els, key=lambda e: e != (1, 0, 0, 1))
    def mul(x, y):
        return ((x[0]*y[0]+x[1]*y[2]) % 3, (x[0]*y[1]+x[1]*y[3]) % 3,
                (x[2]*y[0]+x[3]*y[2]) % 3, (x[2]*y[1]+x[3]*y[3]) % 3)
    return np.array([[idx[mul(a, b)] for b in els] for a in els], dtype=int)


def _gen_dihedral(A):
    """Dih(A) = A : Z2 with the inversion action."""
    m = len(A)
    iv = inverses(A)
    M = 2 * m
    T = np.empty((M, M), dtype=int)
    for x in range(M):
        for y in range(M):
            ix, sx = x % m, x // m
            iy, sy = y % m, y // m
            i = A[ix, iy] if sx == 0 else A[ix, iv[iy]]
            T[x, y] = i + ((sx + sy) % 2) * m
    return T


def _central_product_D4Z4():
    """D4 o Z4, the extraspecial-like group of order 16 (Pauli group)."""
    D = _dihedral(4)
    G = direct(D, cyclic(4))
    # quotient by the diagonal central involution  (z_D, 2)
    o = orders(G)
    zD = None
    for x in range(8):
        if o[x] == 2 and all(D[x, y] == D[y, x] for y in range(8)):
            zD = x
            break
    t = zD * 4 + 2
    K = closure(G, [t])
    cos, rep = {}, []
    for x in range(len(G)):
        c = frozenset(int(G[x, k]) for k in K)
        if c not in cos:
            cos[c] = len(rep)
            rep.append(x)
    M = len(rep)
    lab = {}
    for x in range(len(G)):
        lab[x] = cos[frozenset(int(G[x, k]) for k in K)]
    return np.array([[lab[int(G[rep[a], rep[b]])] for b in range(M)]
                     for a in range(M)], dtype=int)


def _catalogue(limit=34):
    """Named standard groups, keyed by fingerprint."""
    cat = {}

    def put(T, name):
        cat.setdefault(fingerprint2(T), name)
    put(_Sn(3), "S3")
    if limit >= 24:
        put(_Sn(4), "S4")
        put(_SL23(), "SL(2,3)")
    if limit >= 16:
        put(_central_product_D4Z4(), "D4oZ4")
    for n in range(3, limit // 2 + 1):
        put(_dihedral(n), f"D{n}")
    for n in range(2, limit // 4 + 1):
        put(_dicyclic(n), f"Q{4*n}" if (4 * n & (4 * n - 1)) == 0 else f"Dic{n}")
    for k in range(4, 6):
        if 2 ** k <= limit:
            put(_semidihedral(k), f"SD{2**k}")
            put(_modular(k), f"M{2**k}")
    if limit >= 12:
        put(_A4(), "A4")
    # generalised dihedral of every abelian group
    ab = []
    def _ab_tables(n):
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
        for pp in sorted({tuple(sorted(x)) for x in parts(n)}):
            T = None
            for c in pp:
                T = cyclic(c) if T is None else direct(T, cyclic(c))
            yield pp, T
    for n in range(2, limit // 2 + 1):
        for pp, T in _ab_tables(n):
            if len(pp) > 1:
                put(_gen_dihedral(T), "Dih(" + _fmt_abelian(pp) + ")")
    for a in range(3, limit + 1):
        for b in range(2, limit // a + 1):
            if a * b > limit:
                break
            for t in range(2, a):
                if pow(t, b, a) != 1 % a:
                    continue
                T = _semidirect_cyclic(a, b, t)
                put(T, f"Z{a}:Z{b}")
    # direct products of everything named so far with cyclic groups
    for _ in range(2):
        for f, nm in list(cat.items()):
            pass
    return cat


_CAT = None


def name(T):
    global _CAT
    m = len(T)
    if is_abelian(T):
        return _fmt_abelian(_abelian_type(T))
    if _CAT is None:
        _CAT = _catalogue(40)
    f = fingerprint2(T)
    if f in _CAT:
        return _CAT[f]
    for _, n2 in _direct_splits(T):
        return n2
    return None


def _direct_splits(T):
    """G = A x B with A central cyclic; yields (table, name)."""
    m = len(T)
    o = orders(T)
    inv = inverses(T)
    centre = [x for x in range(m)
              if all(T[x, y] == T[y, x] for y in range(m))]
    for z in centre:
        k = int(o[z])
        if k == 1 or k == m:
            continue
        A = closure(T, [z])
        if len(A) != k:
            continue
        # look for a complement that is normal
        for B in _subgroups_of_order(T, m // k):
            if len(A & B) != 1:
                continue
            if not all(T[g, T[b, inv[g]]] in B for g in range(m) for b in B):
                continue
            nb = name(sub_table(T, sorted(B)))
            if nb:
                yield None, f"Z{k}x{nb}"
                return


def _subgroups_of_order(T, k):
    """Subgroups of order k, found from small generating sets."""
    m = len(T)
    out, seen = [], set()
    for a in range(m):
        S = closure(T, [a])
        if len(S) == k and tuple(sorted(S)) not in seen:
            seen.add(tuple(sorted(S)))
            out.append(S)
    for a in range(m):
        for b in range(a + 1, m):
            S = closure(T, [a, b])
            if len(S) == k and tuple(sorted(S)) not in seen:
                seen.add(tuple(sorted(S)))
                out.append(S)
    for a in range(m):
        for b in range(a + 1, m):
            for c in range(b + 1, m):
                S = closure(T, [a, b, c])
                if len(S) == k and tuple(sorted(S)) not in seen:
                    seen.add(tuple(sorted(S)))
                    out.append(S)
    return out


def _heisenberg(p):
    """Extraspecial p^{1+2} of exponent p: upper unitriangular 3x3 over F_p."""
    els = [(a, b, c) for a in range(p) for b in range(p) for c in range(p)]
    els.sort(key=lambda e: e != (0, 0, 0))
    idx = {e: i for i, e in enumerate(els)}
    def mul(x, y):
        return ((x[0] + y[0]) % p, (x[1] + y[1]) % p,
                (x[2] + y[2] + x[0] * y[1]) % p)
    return np.array([[idx[mul(a, b)] for b in els] for a in els], dtype=int)


def names_for(lib):
    """Name every group in the library, with a stable fallback G(n,i)."""
    global _CAT
    if _CAT is None:
        _CAT = _catalogue(40)
    for p in (3, 5):
        _CAT.setdefault(fingerprint2(_heisenberg(p)), f"He{p}")
    out = {}
    for n in sorted(lib):
        used = {}
        for i, T in enumerate(lib[n]):
            nm = name(T) or f"G({n},{i+1})"
            if nm in used:
                used[nm] += 1
                nm = f"{nm}#{used[nm]}"
            else:
                used[nm] = 1
            out[(n, i)] = nm
    return out
