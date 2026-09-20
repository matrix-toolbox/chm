import os, re, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gh_names as NM
import gh_groups as GG

LIB = GG.load_or_build(34)          # cached beside this file as gh_lib.json
NAMES = NM.names_for(LIB)

def inv_table(T):
    n = len(T); e = next(i for i in range(n) if all(T[i][j]==j for j in range(n)))
    inv = [next(j for j in range(n) if T[i][j]==e) for i in range(n)]
    return inv

def quotient(T):
    """Q[i][j] = index of g_j * g_i^{-1}  (the table Appendix F uses)"""
    n = len(T); inv = inv_table(T)
    return np.array([[T[j][inv[i]] for j in range(n)] for i in range(n)])

def product(T):
    return np.array(T)

def load(p, as_phases=True):
    """Octave text blocks.  A real block is read as phases by default -- the
    Catalog stores cores that way -- but with as_phases=False it is taken at
    face value, which is how a real Hadamard matrix of +-1 is stored."""
    txt = open(p).read()
    def blk(name):
        k = txt.find(f'# name: {name}\n')
        if k < 0: return None
        b = txt[k:]
        r = re.search(r'# rows: (\d+)', b)
        if not r: return None
        R, C = int(r.group(1)), int(re.search(r'# columns: (\d+)', b).group(1))
        body = b[re.search(r'# columns: \d+\n', b).end():]
        if 'complex' in b[:b.index('# rows')]:
            t = re.findall(r'\(\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*\)', body)[:R*C]
            return None if len(t)<R*C else np.array([float(x)+1j*float(y) for x,y in t]).reshape(R,C)
        v = re.findall(r'[-\d.eE+]+', body)[:R*C]
        if len(v) < R*C: return None
        A = np.array([float(x) for x in v]).reshape(R, C)
        return A.astype(complex) if not as_phases else np.exp(2j*np.pi*A)
    for nm in ('H','A','B'):
        M = blk(nm)
        if M is not None: return M

def _is_chm(H, tol=1e-6):
    n = H.shape[0]
    return (abs(abs(H) - 1).max() < tol and
            abs(H.conj().T @ H / n - np.eye(n)).max() < tol)


def load_plain(p):
    """A bare array of phases in [0,1), as get_chm writes it: either the full
    N x N array or the (N-1) x (N-1) core.  Whichever is a CHM wins."""
    A = np.loadtxt(p)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        return None
    full = np.exp(2j * np.pi * A)
    if _is_chm(full):
        return full
    n = A.shape[0] + 1                       # treat it as the core instead
    H = np.ones((n, n), dtype=complex)
    H[1:, 1:] = np.exp(2j * np.pi * A)
    return H if _is_chm(H) else None


def load_any(p):
    """.data / .m from the Catalog, or a .gen generator from this appendix."""
    if p.endswith('.gen'):
        import gh_expand
        A = np.asarray(gh_expand.expand(p))
        return A if np.iscomplexobj(A) else np.exp(2j * np.pi * A)
    if '# name:' not in open(p).read(400):
        return load_plain(p)
    H = load(p)
    if H is not None and _is_chm(H):
        return H
    # A real block may already BE the matrix -- a real Hadamard matrix has
    # entries +-1, which load() would otherwise exponentiate as phases.
    raw = load(p, as_phases=False)
    if raw is not None and _is_chm(raw):
        return raw
    # Only accept a block that really is a CHM: a file can carry several blocks,
    # and one read under the wrong interpretation is worse than no answer.
    return None

def dephase(H):
    D = np.diag(1/H[:,0]) @ H
    return D @ np.diag(1/D[0,:])

def developed(M, I):
    n = M.shape[0]; f = np.full(n, np.nan, dtype=complex)
    for i in range(n):
        for j in range(n):
            k = I[i,j]
            if np.isnan(f[k].real): f[k] = M[i,j]
            elif abs(f[k]-M[i,j]) > 1e-6: return False
    return True

def same_multisets(M, tol=1e-5):
    """In M=[f(g_i g_j)] each row is f applied to a permutation of G, so every
    row carries one value multiset, and likewise every column.  Relabelling the
    group only permutes rows and columns, so failing this proves M is not
    group-developed over ANY group, in any ordering -- unlike classify(), whose
    negative answer only concerns the stored ordering.

    The entries are clustered once by circular distance and compared as
    multisets of cluster indices.  Rounding the phase instead puts values either
    side of the 0 / 2*pi wrap into different bins, so two entries a femtoradian
    apart would be called different.
    """
    n = M.shape[0]
    reps, idx = [], np.empty(M.size, dtype=int)
    flat = M.ravel()
    for t, z in enumerate(flat):
        for k, r in enumerate(reps):
            if abs(z - r) < tol:
                idx[t] = k
                break
        else:
            reps.append(z)
            idx[t] = len(reps) - 1
    I = idx.reshape(M.shape)
    rows = {tuple(sorted(I[i])) for i in range(n)}
    cols = {tuple(sorted(I[:, j])) for j in range(n)}
    return len(rows) == 1 and len(cols) == 1


def recover_group(C, tol=1e-6):
    """Is the core C[g,h] = F(h g^-1) for SOME group law, under any labelling?

    classify() tests one fixed ordering of the group elements and its negative
    answer therefore means nothing on its own -- a matrix developed over G but
    stored in a different order comes back "neither".  This does not depend on
    the ordering: fix row 0, and for every row find the column permutation
    carrying it onto row 0.  Those are the regular representation of G, so they
    must be distinct and closed under composition.  Needs F injective; when it
    is not, the answer is 'inconclusive', never 'no'.
    """
    n = C.shape[0]
    perms = []
    for i in range(n):
        p = [-1] * n
        for j in range(n):
            hit = [k for k in range(n) if abs(C[i, j] - C[0, k]) < tol]
            if len(hit) != 1:
                return None, 'inconclusive: generator not injective'
            p[j] = hit[0]
        if sorted(p) != list(range(n)):
            return None, 'no: a row is not a permutation of row 0'
        perms.append(tuple(p))
    if len(set(perms)) != n:
        return None, 'no: permutations not distinct'
    S = set(perms)
    for a in perms:
        for b in perms:
            if tuple(a[b[k]] for k in range(n)) not in S:
                return None, 'no: not closed under composition'
    ab = all(tuple(a[b[k]] for k in range(n)) == tuple(b[a[k]] for k in range(n))
             for a in perms for b in perms)
    return perms, 'yes: %s group of order %d' % ('abelian' if ab else 'non-abelian', n)


def classify(H):
    N = H.shape[0]
    core = dephase(H)[1:,1:]
    for gi, T in enumerate(LIB.get(N-1, [])):
        for tag, I in (('quotient', quotient(T)), ('product', product(T))):
            for M, t2 in ((core,''), (core.T,' transposed')):
                if developed(M, I): return f'BORDERED {tag} / {NAMES[(N-1,gi)]}{t2}'
    for gi, T in enumerate(LIB.get(N, [])):
        for tag, I in (('quotient', quotient(T)), ('product', product(T))):
            for M, t2 in ((H,''), (dephase(H),' dephased')):
                if developed(M, I): return f'UNBORDERED {tag} / {NAMES[(N,gi)]}{t2}'
    return 'not found in the tested orderings'


if __name__ == "__main__":
    # usage: python3 gh_check.py FILE...   -- is each matrix group-developed?
    #
    # Two answers are printed, because they scope differently.  classify() tries
    # the cached group tables in their stored element order, so its negative is
    # only "not in the orderings tested" -- and the cache stops at order 34.
    # recover_group() reads the group off the matrix itself and does not depend
    # on any ordering; it says "no" only on a real structural failure.
    import glob
    for pat in sys.argv[1:]:
        for p in sorted(glob.glob(pat)) or [pat]:
            try:
                H = load_any(p)
            except Exception as e:
                print(f"{p:<44} unreadable: {e}")
                continue
            if H is None:
                print(f"{p:<44} unreadable")
                continue
            D = H / (H[:, :1] * H[:1, :] / H[0, 0])
            _, why = recover_group(D[1:, 1:])
            print(f"{p:<44} {classify(H)}")
            print(f"{'':<44}   ordering-free: {why}")
