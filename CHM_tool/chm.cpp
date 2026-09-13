// =============================================================================
// chm.cpp -- Complex Hadamard Matrices: shared core library (implementation)
//
// 2026-09-04  W. Bruzda  <w.bruzda@cft.edu.pl>
// See chm.hpp for the interface and COMPILE.txt for build lines.
// =============================================================================

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L      // popen/pclose are POSIX, not ISO C++
#endif

// Make LAPACKE speak std::complex<double> directly, so no reinterpret_cast or
// buffer copying is ever needed.  lapack.h honours a pre-defined
// lapack_complex_double; LAPACK_COMPLEX_CUSTOM stops lapacke_config.h from
// redefining it to "double _Complex" afterwards.  (This is the recipe
// documented in the comment at the top of /usr/include/lapack.h.)
#include <complex>
#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float  std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>

#define CHM_PROGRAM_ABI 20260908   // must match CHM_ABI in chm.hpp
#include "chm.hpp"
#if CHM_ABI < CHM_PROGRAM_ABI
#error "MIXED SOURCE TREE: chm.hpp is older than this .cpp. Delete the directory and extract one archive into it; see COMPILE.txt."
#endif


#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <unordered_map>
#include <cstdlib>           // strtoll
#include <unistd.h>          // readlink, for locating the sibling programs

#ifdef _OPENMP
#include <omp.h>
#endif

// OpenBLAS lets us pin its internal thread pool.  Declared weak so that the
// code still links against the reference LAPACK/BLAS, where it is absent.
extern "C" void openblas_set_num_threads(int) __attribute__((weak));

namespace chm {

namespace {

inline void blas_threads(int t) {
    if (openblas_set_num_threads) openblas_set_num_threads(t);
}

// Row-major C = A * B, hand-rolled: N is small (typically 4..64) and this keeps
// the library free of a CBLAS header dependency.  The i-k-j order streams both
// B and C linearly, which vectorises cleanly under -O3 -march=znver3.
CMat mul(const CMat& A, const CMat& B) {
    const int n = A.n;
    CMat C(n);
    for (int i = 0; i < n; ++i) {
        const cd* arow = &A.a[std::size_t(i) * n];
        cd*       crow = &C.a[std::size_t(i) * n];
        for (int k = 0; k < n; ++k) {
            const cd  aik = arow[k];
            if (aik == cd(0.0, 0.0)) continue;
            const cd* brow = &B.a[std::size_t(k) * n];
            for (int j = 0; j < n; ++j) crow[j] += aik * brow[j];
        }
    }
    return C;
}

// Row-major Gram matrix G = X * X^dagger (Hermitian; only |G| matters to us).
CMat gram(const CMat& X) {
    const int n = X.n;
    CMat G(n);
    for (int i = 0; i < n; ++i) {
        const cd* xi = &X.a[std::size_t(i) * n];
        for (int j = i; j < n; ++j) {
            const cd* xj = &X.a[std::size_t(j) * n];
            cd s(0.0, 0.0);
            for (int k = 0; k < n; ++k) s += xi[k] * std::conj(xj[k]);
            G(i, j) = s;
            if (i != j) G(j, i) = std::conj(s);
        }
    }
    return G;
}

inline double wrap01(double x) {
    x = std::fmod(x, 1.0);
    if (x < 0.0) x += 1.0;
    return x;
}

// arg(z) mapped to [0,1)
inline double phase01(cd z) { return wrap01(std::arg(z) / TAU); }

inline bool is_ident_char(char c, bool first) {
    if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') return true;
    return !first && std::isdigit(static_cast<unsigned char>(c));
}

}  // namespace

// -----------------------------------------------------------------------------
// CMat helpers
// -----------------------------------------------------------------------------
CMat CMat::identity(int n) {
    CMat I(n);
    for (int i = 0; i < n; ++i) I(i, i) = cd(1.0, 0.0);
    return I;
}

CMat CMat::fourier(int n) {
    CMat F(n);
    for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k)
            F(j, k) = std::polar(1.0, TAU * (double(j) * k) / n);
    return F;
}

// -----------------------------------------------------------------------------
// Norms and elementary transformations
// -----------------------------------------------------------------------------
double n1(const CMat& X) {
    double s = 0.0;
    for (const cd& z : X.a) { const double d = std::abs(z) - 1.0; s += d * d; }
    return std::sqrt(s);
}

double nh(const CMat& X) {
    const int n = X.n;
    const CMat G = gram(X);
    double s = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            cd d = G(i, j);
            if (i == j) d -= cd(double(n), 0.0);
            s += std::norm(d);
        }
    return std::sqrt(s);
}

double frob(const CMat& X) {
    double s = 0.0;
    for (const cd& z : X.a) s += std::norm(z);
    return std::sqrt(s);
}

CMat transpose(const CMat& X) {
    CMat Y(X.n);
    for (int i = 0; i < X.n; ++i)
        for (int j = 0; j < X.n; ++j) Y(i, j) = X(j, i);
    return Y;
}

CMat adjoint(const CMat& X) {
    CMat Y(X.n);
    for (int i = 0; i < X.n; ++i)
        for (int j = 0; j < X.n; ++j) Y(i, j) = std::conj(X(j, i));
    return Y;
}

double sym_error(const CMat& X) {
    double s = 0.0;
    for (int i = 0; i < X.n; ++i)
        for (int j = 0; j < X.n; ++j) s += std::norm(X(i, j) - X(j, i));
    return std::sqrt(s);
}

double herm_error(const CMat& X) {
    double s = 0.0;
    for (int i = 0; i < X.n; ++i)
        for (int j = 0; j < X.n; ++j) s += std::norm(X(i, j) - std::conj(X(j, i)));
    return std::sqrt(s);
}

// Port of dephase.m: K = diag(1/H(:,1)) * H * diag(1/H(1,:)), with the border
// forced to exact ones.  This convention preserves both H = H^T and H = H^dag.
CMat dephase(const CMat& H) {
    const int n = H.n;
    CMat K(n);
    for (int j = 0; j < n; ++j) {
        const cd cj = H(0, j);
        for (int i = 0; i < n; ++i) K(i, j) = H(i, j) / cj;
    }
    for (int i = 0; i < n; ++i) {
        const cd ri = K(i, 0);
        for (int j = 0; j < n; ++j) K(i, j) /= ri;
    }
    for (int j = 0; j < n; ++j) K(0, j) = cd(1.0, 0.0);
    for (int i = 0; i < n; ++i) K(i, 0) = cd(1.0, 0.0);
    return K;
}

// Port of reshuffle.m:  A^R_{(m,a),(n,b)} = A_{(m,n),(a,b)}   with N = d*d.
CMat reshuffle(const CMat& X) {
    const int n = X.n;
    const int d = int(std::lround(std::sqrt(double(n))));
    if (d * d != n) return CMat();
    CMat Y(n);
    for (int m = 0; m < d; ++m)
        for (int a = 0; a < d; ++a)
            for (int nn = 0; nn < d; ++nn)
                for (int b = 0; b < d; ++b)
                    Y(m * d + a, nn * d + b) = X(m * d + nn, a * d + b);
    return Y;
}

// Port of Tx.m for a bipartite d x d split.
//   sys = 2 :  A^{G}_{(m,a),(n,b)} = A_{(m,b),(n,a)}
//   sys = 1 :  A^{G}_{(m,a),(n,b)} = A_{(n,a),(m,b)}
// The two differ by a global transpose, so every unitarily invariant quantity
// (in particular the linear entropy) agrees on them.
CMat partial_transpose(const CMat& X, int sys) {
    const int n = X.n;
    const int d = int(std::lround(std::sqrt(double(n))));
    if (d * d != n) return CMat();
    CMat Y(n);
    for (int m = 0; m < d; ++m)
        for (int a = 0; a < d; ++a)
            for (int nn = 0; nn < d; ++nn)
                for (int b = 0; b < d; ++b) {
                    Y(m * d + a, nn * d + b) =
                        (sys == 2) ? X(m * d + b, nn * d + a)
                                   : X(nn * d + a, m * d + b);
                }
    return Y;
}

std::vector<std::vector<double>> core_phases(const CMat& H) {
    const int m = H.n - 1;
    std::vector<std::vector<double>> A(std::max(m, 0), std::vector<double>(std::max(m, 0), 0.0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j) A[i][j] = phase01(H(i + 1, j + 1));
    return A;
}

CMat from_core_phases(const std::vector<std::vector<double>>& A) {
    const int m = int(A.size());
    const int n = m + 1;
    CMat H(n, cd(1.0, 0.0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j) H(i + 1, j + 1) = std::polar(1.0, TAU * A[i][j]);
    return H;
}

CMat from_full_phases(const std::vector<std::vector<double>>& A) {
    const int n = int(A.size());
    CMat H(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) H(i, j) = std::polar(1.0, TAU * A[i][j]);
    return H;
}

// -----------------------------------------------------------------------------
// Linear algebra: LAPACKE backend
// -----------------------------------------------------------------------------
// Singular values only.  A row-major array of shape (rows x cols) and leading
// dimension cols is, read as column-major with the dimensions swapped, exactly
// A^T -- and A and A^T share their singular values.  Handing LAPACKE the
// column-major view therefore costs nothing and, unlike LAPACK_ROW_MAJOR,
// needs no U/V^dag workspace at all for jobz = 'N' (the row-major wrapper
// rejects ldvt = 1 there).
std::vector<double> singular_values(const RMat& A) {
    const lapack_int m = A.cols, n = A.rows;        // dimensions of A^T
    const lapack_int k = std::min(m, n);
    std::vector<double> a = A.a;                    // gesdd overwrites its input
    std::vector<double> s(std::size_t(k), 0.0);
    double du = 0.0, dvt = 0.0;                     // unreferenced for jobz = 'N'
    const lapack_int info = LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'N', m, n, a.data(), m,
                                           s.data(), &du, 1, &dvt, 1);
    if (info != 0) {
        std::fprintf(stderr, "chm: LAPACKE_dgesdd failed, info = %d\n", int(info));
        s.assign(std::size_t(k), std::nan(""));
    }
    return s;
}

std::vector<double> singular_values(const CMat& A) {
    const lapack_int n = A.n;
    std::vector<cd> a = A.a;
    std::vector<double> s(std::size_t(n), 0.0);
    cd du(0.0, 0.0), dvt(0.0, 0.0);
    const lapack_int info = LAPACKE_zgesdd(LAPACK_COL_MAJOR, 'N', n, n, a.data(), n,
                                           s.data(), &du, 1, &dvt, 1);
    if (info != 0) {
        std::fprintf(stderr, "chm: LAPACKE_zgesdd failed, info = %d\n", int(info));
        s.assign(std::size_t(n), std::nan(""));
    }
    return s;
}

// Unitary polar factor: X = W S V^dag  ->  U = W V^dag.
// Replaces PD.m / PD_SVD.m / PD_SVDi.m.  The SVD route is used rather than the
// eig(X^dag X) route of PD_SVDi.m because the latter squares the condition
// number and loses about half the significant digits near singular X.
CMat polar_unitary(const CMat& X) {
    const lapack_int n = X.n;
    std::vector<cd> a = X.a;
    std::vector<double> s(std::size_t(n), 0.0);
    CMat W{int(n)}, Vt{int(n)};
    const lapack_int info = LAPACKE_zgesdd(LAPACK_ROW_MAJOR, 'S', n, n, a.data(), n, s.data(),
                                           W.a.data(), n, Vt.a.data(), n);
    if (info != 0) return CMat();
    return mul(W, Vt);
}

int rank_from_sv(const std::vector<double>& s, int rows, int cols, double tol) {
    if (s.empty()) return 0;
    double smax = 0.0;
    for (double v : s) smax = std::max(smax, v);
    const double cut = (tol > 0.0)
                           ? tol
                           : double(std::max(rows, cols)) * 2.220446049250313e-16 * smax;
    int r = 0;
    for (double v : s) if (v > cut) ++r;
    return r;
}

// -----------------------------------------------------------------------------
// Defect  (port of ud.m)
// -----------------------------------------------------------------------------
namespace {

// System matrix R of the linearised phase-perturbation conditions.
// Rows: 2 per pair (row < next_row), split into real and imaginary part.
// Columns: the N*N unknown phase increments.
RMat build_R(const CMat& U) {
    const int n = U.n;
    const int tau = n * (n - 1) / 2;
    RMat R(2 * tau, n * n);
    int p = 0;
    for (int r = 0; r < n - 1; ++r) {
        for (int q = r + 1; q < n; ++q, ++p) {
            const int t0 = 2 * p, t1 = 2 * p + 1;
            for (int k = 0; k < n; ++k) {
                const cd M = U(r, k) * std::conj(U(q, k));
                R(t0, r * n + k) = M.real();
                R(t1, r * n + k) = M.imag();
                R(t0, q * n + k) = -M.real();
                R(t1, q * n + k) = -M.imag();
            }
        }
    }
    return R;
}

}  // namespace

int defect(const CMat& U, char method, double tol) {
    const int n = U.n;
    if (n < 1) return 0;
    // Method 'T' of ud.m builds the matrix of images of tangent vectors, which
    // is exactly the transpose of R above (entry by entry).  Since a matrix and
    // its transpose share their rank, 'T' and 'R' are mathematically identical
    // and are handled by the same code path.
    const RMat R = build_R(U);
    const std::vector<double> s = singular_values(R);
    const double use_tol = (method == 'S') ? tol : -1.0;
    const int rk = rank_from_sv(s, R.rows, R.cols, use_tol);
    return (n - 1) * (n - 1) - rk;
}

// -----------------------------------------------------------------------------
// Haagerup invariants  (replaces cL / getUnique.m)
// -----------------------------------------------------------------------------
namespace {

// Tolerance-based uniqueness via a uniform hash grid of cell size eps.
// getUnique.m does the same job in O(M^2); the grid makes it O(M) expected,
// which matters because M = O(N^4) candidates are generated.
class UniqueSet {
  public:
    explicit UniqueSet(double eps) : eps_(eps > 0 ? eps : 1e-12) {}

    bool insert(cd z) {
        if (!reps_.empty() && std::abs(z - reps_[last_]) < eps_) return false;  // locality
        const long long cx = llround(std::floor(z.real() / eps_));
        const long long cy = llround(std::floor(z.imag() / eps_));
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy) {
                auto it = grid_.find(key(cx + dx, cy + dy));
                if (it == grid_.end()) continue;
                for (int idx : it->second)
                    if (std::abs(z - reps_[std::size_t(idx)]) < eps_) { last_ = std::size_t(idx); return false; }
            }
        reps_.push_back(z);
        last_ = reps_.size() - 1;
        grid_[key(cx, cy)].push_back(int(reps_.size()) - 1);
        return true;
    }

    const std::vector<cd>& reps() const { return reps_; }

  private:
    static std::uint64_t key(long long x, long long y) {
        return (std::uint64_t(std::uint32_t(x)) << 32) | std::uint64_t(std::uint32_t(y));
    }
    double eps_;
    std::unordered_map<std::uint64_t, std::vector<int>> grid_;
    std::vector<cd> reps_;
    std::size_t last_ = 0;
};

}  // namespace

// Lambda(H) = { H_ij H_kl / (H_il H_kj) }.  With r_j = H_ij / H_kj the quadruple
// invariant factorises as r_j / r_l, which cuts the work from O(N^4) divisions
// with four loads each to O(N^2) ratios per row pair.  Row pairs (i,k) and (k,i)
// produce reciprocal sets and (j,l)/(l,j) already covers reciprocals, so only
// i < k has to be visited.
std::vector<cd> haagerup_set(const CMat& H, double eps) {
    const int n = H.n;
    std::vector<cd> out;
    if (n < 2) { out.push_back(cd(1.0, 0.0)); return out; }

    struct Pair { int i, k; };
    std::vector<Pair> pairs;
    pairs.reserve(std::size_t(n) * (n - 1) / 2);
    for (int i = 0; i < n - 1; ++i)
        for (int k = i + 1; k < n; ++k) pairs.push_back({i, k});

    std::vector<std::vector<cd>> partial;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        UniqueSet local(eps);
        std::vector<cd> r((std::size_t(n)));
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (std::size_t pi = 0; pi < pairs.size(); ++pi) {
            const int i = pairs[pi].i, k = pairs[pi].k;
            for (int j = 0; j < n; ++j) r[std::size_t(j)] = H(i, j) / H(k, j);
            for (int j = 0; j < n; ++j)
                for (int l = 0; l < n; ++l) local.insert(r[std::size_t(j)] / r[std::size_t(l)]);
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        partial.push_back(local.reps());
    }

    UniqueSet merged(eps);
    merged.insert(cd(1.0, 0.0));  // i = k and j = l always contribute unity
    for (const auto& v : partial)
        for (cd z : v) merged.insert(z);

    out = merged.reps();
    std::sort(out.begin(), out.end(), [](cd a, cd b) {
        const double aa = std::arg(a), ab = std::arg(b);
        if (aa != ab) return aa < ab;
        return std::abs(a) < std::abs(b);
    });
    return out;
}

// -----------------------------------------------------------------------------
// Butson class  (port of isBH.m and LF.m)
// -----------------------------------------------------------------------------
// isBH.m raises the entries to the power q, which accumulates round-off for
// large q.  Testing the phases instead is exact to within one rounding:
//   H_jk^q = 1  <=>  q * theta_jk in Z,  and  |exp(2 pi i q theta) - 1| = 2|sin(pi q theta)|.
int butson_q(const CMat& H, int qmax, double tol) {
    const int n = H.n;
    std::vector<double> th(H.a.size());
    for (std::size_t i = 0; i < H.a.size(); ++i) th[i] = phase01(H.a[i]);
    for (int q = 2; q <= qmax; ++q) {
        double worst = 0.0;
        for (std::size_t i = 0; i < th.size(); ++i) {
            const double x = double(q) * th[i];
            const double e = 2.0 * std::abs(std::sin(PI * (x - std::floor(x + 0.5))));
            if (e > worst) { worst = e; if (worst > tol) break; }
        }
        if (worst <= tol) return q;
    }
    (void)n;
    return 0;
}

std::vector<std::vector<int>> log_form(const CMat& H, int q) {
    const int n = H.n;
    std::vector<std::vector<int>> L(std::size_t(n), std::vector<int>(std::size_t(n), 0));
    if (q <= 0) return L;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            long long v = llround(phase01(H(i, j)) * q) % q;
            if (v < 0) v += q;
            L[std::size_t(i)][std::size_t(j)] = int(v);
        }
    return L;
}

// -----------------------------------------------------------------------------
// Linear entropy  (port of SL.m and SL3.m)
// -----------------------------------------------------------------------------
// f = M M^dag / Tr(M M^dag);  SL = (1 - Tr f^2) * N/(N-1);  SL(unitary) = 1.
// f is Hermitian, so Tr f^2 = sum |f_ab|^2 -- no eigendecomposition needed.
double linear_entropy(const CMat& M) {
    const int n = M.n;
    if (n < 2) return 0.0;
    const CMat G = gram(M);
    double tr = 0.0;
    for (int i = 0; i < n; ++i) tr += G(i, i).real();
    if (tr <= 0.0) return 0.0;
    double s2 = 0.0;
    for (const cd& z : G.a) s2 += std::norm(z);
    const double trf2 = s2 / (tr * tr);
    return (1.0 - trf2) * double(n) / double(n - 1);
}

SL3 linear_entropy_triplet(const CMat& H, double tol) {
    SL3 r;
    const int n = H.n;
    const int d = int(std::lround(std::sqrt(double(n))));
    if (d * d != n || d < 2) return r;  // H^R and H^G need a bipartite d x d split
    r.defined = true;
    r.d       = d;
    r.s_H     = linear_entropy(H);
    r.s_R     = linear_entropy(reshuffle(H));
    r.s_G     = linear_entropy(partial_transpose(H, 2));
    r.r_dual  = std::abs(r.s_H - 1.0) < tol && std::abs(r.s_R - 1.0) < tol;
    r.g_dual  = std::abs(r.s_H - 1.0) < tol && std::abs(r.s_G - 1.0) < tol;
    r.two_unitary = r.r_dual && r.g_dual;
    return r;
}

// -----------------------------------------------------------------------------
// Generation
// -----------------------------------------------------------------------------
namespace {

inline void project_symmetry(CMat& X, Symmetry sym) {
    const int n = X.n;
    if (sym == Symmetry::Symmetric) {
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                const cd m = 0.5 * (X(i, j) + X(j, i));
                X(i, j) = m; X(j, i) = m;
            }
    } else if (sym == Symmetry::Hermitian) {
        for (int i = 0; i < n; ++i) {
            X(i, i) = cd(X(i, i).real(), 0.0);
            for (int j = i + 1; j < n; ++j) {
                const cd m = 0.5 * (X(i, j) + std::conj(X(j, i)));
                X(i, j) = m; X(j, i) = std::conj(m);
            }
        }
    }
}

inline void unimodularise(CMat& X) {
    for (cd& z : X.a) {
        const double m = std::abs(z);
        z = (m > 0.0) ? z / m : cd(1.0, 0.0);
    }
}

// --- Sinkhorn (port of sinkhorn.m / ss.m / sh.m) ------------------------------
// Alternating projection between the unimodular matrices and the sphere of
// radius sqrt(N) unitaries.  The map has fixed points that are not CHM, and for
// N beyond roughly 12 a plain run stagnates at n1 ~ 1e-2 (exactly the behaviour
// the "Z" bookkeeping in sinkhorn.m/ss.m tries to detect).  So when progress
// stops we first perturb the phases a few times -- much cheaper than throwing
// the iterate away -- and only then report the restart as failed.
bool sinkhorn_restart(CMat& X, const GenOptions& opt, std::mt19937_64& rng,
                      long long& iters_used, double& v1, double& vh) {
    const int n = opt.size;
    std::uniform_real_distribution<double> U(0.0, 1.0);

    X = CMat(n);
    for (cd& z : X.a) z = std::polar(1.0, TAU * U(rng));
    project_symmetry(X, opt.symmetry);

    const double sq = std::sqrt(double(n));
    double    best  = std::numeric_limits<double>::infinity();
    long long stall = 0;
    int       kicks_left = std::max(0, opt.kicks);
    const long long stall_limit = std::max<long long>(64, opt.iterations / 8);
    v1 = vh = std::numeric_limits<double>::infinity();

    for (iters_used = 0; iters_used < opt.iterations; ++iters_used) {
        if (v1 < opt.eps && vh < opt.eps) return true;

        unimodularise(X);                       // enforce |H_jk| = 1
        CMat P = polar_unitary(X);              // nearest unitary
        if (P.n == 0) return false;             // LAPACK failure
        for (cd& z : P.a) z *= sq;              // scale to H H^dag = N I
        X = std::move(P);
        project_symmetry(X, opt.symmetry);      // re-impose H = H^T or H = H^dag

        v1 = n1(X);
        vh = nh(X);

        if (v1 < best * (1.0 - 1e-14)) {
            best  = v1;
            stall = 0;
        } else if (++stall > stall_limit) {
            if (kicks_left-- <= 0) break;
            // Random phase kick, decreasing in size, keeping the symmetry.
            const double sigma = 0.15 / double(1 + opt.kicks - kicks_left);
            for (cd& z : X.a) z *= std::polar(1.0, TAU * sigma * (2.0 * U(rng) - 1.0));
            project_symmetry(X, opt.symmetry);
            best  = std::numeric_limits<double>::infinity();
            stall = 0;
        }
    }
    return (v1 < opt.eps && vh < opt.eps);
}

// --- rwcp : random walk over core phases --------------------------------------
// State is the (N-1)x(N-1) array of core phases A in [0,1); the matrix
//     H = [ 1 ... 1 ; 1 exp(2 pi i A) ]
// is unimodular and dephased by construction, so n1(H) = 0 exactly and the only
// objective left is
//     F(A) = nh(H)^2 = sum_{a != b} |G_ab|^2 ,   G = H H^dag
// (the diagonal is exactly N because every entry has modulus one).
//
// Moving a single phase changes exactly one entry of H, hence only one row and
// one column of G: the objective can be updated in O(N) instead of O(N^3).
//
// Two move types are drawn at random:
//   * an exploratory step  A_pq += uniform(-step, step)  with an adaptive step;
//   * the exact coordinate minimiser.  Writing G_rb = C_rb + u * conj(H_bc) with
//     u = H_rc on the unit circle, the r-part of F is
//         const + 4 Re( u * Z ),   Z = sum_{b != r} conj(C_rb) conj(H_bc),
//     which is minimal at u = -Z/|Z| -- available in O(N) as well.
// Every proposal is scored by its exact dF and accepted only if it does not
// increase F, except for the deliberate kicks used to escape a stalled walk.
class RWCP {
  public:
    RWCP(const GenOptions& o, std::mt19937_64& rng)
        : opt_(o), rng_(rng), n_(o.size), m_(o.size - 1) {}

    bool restart(long long& iters_used, double& v1, double& vh) {
        const double target = opt_.eps * opt_.eps;
        std::uniform_real_distribution<double> U(0.0, 1.0);
        const int md = std::max(m_, 1);

        randomise();
        rebuild();

        double    step  = 0.25;
        double    best  = F_;
        double    best_ever = F_;
        long long stall = 0;
        const long long stall_limit = std::max<long long>(256, 64LL * md * md);

        for (iters_used = 0; iters_used < opt_.iterations; ++iters_used) {
            if (F_ <= target) break;

            const int  p     = int(U(rng_) * md) % md;
            const int  q     = int(U(rng_) * md) % md;
            const bool exact = (U(rng_) < 0.75);
            const double cand = exact ? exact_min_phase(p + 1, q + 1)
                                      : wrap01(A_[idx(p, q)] + step * (2.0 * U(rng_) - 1.0));

            if (opt_.symmetry == Symmetry::None) try_free(p, q, cand, false);
            else                                 try_mirrored(p, q, cand, false);

            if (F_ < best_ever) best_ever = F_;
            if (F_ < best * (1.0 - 1e-13)) { best = F_; stall = 0; }
            else if (++stall > stall_limit) {
                // Annealing kick: force a few phases to fresh random values,
                // accepting the uphill move, then continue descending.
                const int kicks = 1 + int(U(rng_) * md);
                for (int t = 0; t < kicks; ++t) {
                    const int kp = int(U(rng_) * md) % md;
                    const int kq = int(U(rng_) * md) % md;
                    if (opt_.symmetry == Symmetry::None) try_free(kp, kq, U(rng_), true);
                    else                                 try_mirrored(kp, kq, U(rng_), true);
                }
                best  = F_;
                stall = 0;
                step  = 0.25;
            }

            // Shrink the exploratory step slowly and refresh G from H now and
            // then, so that incremental round-off cannot accumulate.
            if ((iters_used & 0x3FF) == 0x3FF) {
                step = std::max(1e-17, step * 0.9);
                rebuild();
            }
        }

        rebuild();
        v1 = n1(H_);                                        // exactly 0 by construction
        const bool ok = (v1 < opt_.eps && std::sqrt(std::max(0.0, F_)) < opt_.eps);
        vh = ok ? std::sqrt(std::max(0.0, F_)) : std::sqrt(std::max(0.0, best_ever));
        return ok;
    }

    const CMat& matrix() const { return H_; }

    // Seed the walk from an existing matrix instead of from random phases: the
    // moduli are discarded (that is the unimodular projection), the result is
    // dephased and the symmetry constraint is re-imposed on the core.  Used to
    // polish a Sinkhorn iterate that has come close but converges only linearly.
    void seed_from(const CMat& H) {
        const CMat K = dephase(H);
        A_.assign(std::size_t(std::max(m_, 1)) * std::max(m_, 1), 0.0);
        for (int p = 0; p < m_; ++p)
            for (int q = 0; q < m_; ++q) A_[idx(p, q)] = phase01(K(p + 1, q + 1));

        if (opt_.symmetry == Symmetry::Symmetric) {
            for (int p = 0; p < m_; ++p)
                for (int q = p + 1; q < m_; ++q) {
                    // circular mean of the two phases
                    const cd m = std::polar(1.0, TAU * A_[idx(p, q)]) +
                                 std::polar(1.0, TAU * A_[idx(q, p)]);
                    const double a = (std::abs(m) > 1e-12) ? phase01(m) : A_[idx(p, q)];
                    A_[idx(p, q)] = a;
                    A_[idx(q, p)] = a;
                }
        } else if (opt_.symmetry == Symmetry::Hermitian) {
            for (int p = 0; p < m_; ++p) {
                A_[idx(p, p)] = (A_[idx(p, p)] > 0.25 && A_[idx(p, p)] < 0.75) ? 0.5 : 0.0;
                for (int q = p + 1; q < m_; ++q) {
                    const cd m = std::polar(1.0, TAU * A_[idx(p, q)]) +
                                 std::polar(1.0, TAU * wrap01(-A_[idx(q, p)]));
                    const double a = (std::abs(m) > 1e-12) ? phase01(m) : A_[idx(p, q)];
                    A_[idx(p, q)] = a;
                    A_[idx(q, p)] = wrap01(-a);
                }
            }
        }
        rebuild();
    }

    // Pure exact-coordinate descent, no exploration.  Locally the objective is
    // a well-behaved quadratic, so this converges to machine precision quickly
    // once the iterate is in the right basin.
    bool sharpen(long long iters, double& v1, double& vh) {
        const double target = opt_.eps * opt_.eps;
        std::uniform_real_distribution<double> U(0.0, 1.0);
        const int md = std::max(m_, 1);
        for (long long it = 0; it < iters && F_ > target; ++it) {
            const int p = int(U(rng_) * md) % md;
            const int q = int(U(rng_) * md) % md;
            const double cand = exact_min_phase(p + 1, q + 1);
            if (opt_.symmetry == Symmetry::None) try_free(p, q, cand, false);
            else                                 try_mirrored(p, q, cand, false);
            if ((it & 0xFFF) == 0xFFF) rebuild();
        }
        rebuild();
        v1 = n1(H_);
        vh = std::sqrt(std::max(0.0, F_));
        return (v1 < opt_.eps && vh < opt_.eps);
    }

  private:
    std::size_t idx(int p, int q) const { return std::size_t(p) * std::max(m_, 1) + q; }

    void randomise() {
        std::uniform_real_distribution<double> U(0.0, 1.0);
        A_.assign(std::size_t(std::max(m_, 1)) * std::max(m_, 1), 0.0);
        for (int p = 0; p < m_; ++p)
            for (int q = 0; q < m_; ++q) A_[idx(p, q)] = U(rng_);
        if (opt_.symmetry == Symmetry::Symmetric) {
            for (int p = 0; p < m_; ++p)
                for (int q = p + 1; q < m_; ++q) A_[idx(q, p)] = A_[idx(p, q)];
        } else if (opt_.symmetry == Symmetry::Hermitian) {
            // A dephased Hermitian matrix has  A_qp = -A_pq (mod 1)  off the
            // diagonal and  A_pp in {0, 1/2}  on it (H_pp real and unimodular).
            for (int p = 0; p < m_; ++p) {
                A_[idx(p, p)] = (U(rng_) < 0.5) ? 0.0 : 0.5;
                for (int q = p + 1; q < m_; ++q) A_[idx(q, p)] = wrap01(-A_[idx(p, q)]);
            }
        }
    }

    void rebuild() {
        H_ = CMat(n_, cd(1.0, 0.0));
        for (int p = 0; p < m_; ++p)
            for (int q = 0; q < m_; ++q) H_(p + 1, q + 1) = std::polar(1.0, TAU * A_[idx(p, q)]);
        G_ = gram(H_);
        F_ = objective();
    }

    double objective() const {
        double s = 0.0;
        for (int a = 0; a < n_; ++a)
            for (int b = 0; b < n_; ++b)
                if (a != b) s += std::norm(G_(a, b));
        return s;
    }

    // Recompute rows r of G from H (and mirror into the columns): O(N^2).
    void refresh_row(int r) {
        for (int b = 0; b < n_; ++b) {
            cd s(0.0, 0.0);
            for (int k = 0; k < n_; ++k) s += H_(r, k) * std::conj(H_(b, k));
            G_(r, b) = s;
            G_(b, r) = std::conj(s);
        }
    }

    // ---- unconstrained move: one entry, exact O(N) update --------------------
    bool try_free(int p, int q, double phase, bool allow_uphill) {
        const int r = p + 1, c = q + 1;
        const cd unew = std::polar(1.0, TAU * phase);
        const cd du   = unew - H_(r, c);

        // Only row r (and, by Hermiticity, column r) of G changes; the terms of
        // F touching row r contribute 2 * sum_{b != r} |G_rb|^2.
        double d = 0.0;
        for (int b = 0; b < n_; ++b) {
            if (b == r) continue;
            const cd g = G_(r, b) + du * std::conj(H_(b, c));
            d += 2.0 * (std::norm(g) - std::norm(G_(r, b)));
        }
        if (!allow_uphill && d > 0.0) return false;

        H_(r, c) = unew;
        for (int b = 0; b < n_; ++b) {
            if (b == r) continue;
            const cd g = G_(r, b) + du * std::conj(H_(b, c));
            G_(r, b) = g;
            G_(b, r) = std::conj(g);
        }
        A_[idx(p, q)] = phase;
        F_ = std::max(0.0, F_ + d);
        return true;
    }

    // ---- constrained move: entry plus its mirror, O(N^2) ---------------------
    // Under H = H^T or H = H^dag the partner entry must move in lock step, so
    // the pair is scored exactly and rolled back when it would raise F.
    bool try_mirrored(int p, int q, double phase, bool allow_uphill) {
        const double old_pq = A_[idx(p, q)];
        const double old_qp = A_[idx(q, p)];
        const double f_before = F_;

        set_mirrored(p, q, phase);
        const int r1 = p + 1, r2 = q + 1;
        refresh_row(r1);
        if (r2 != r1) refresh_row(r2);
        F_ = objective();

        if (allow_uphill || F_ <= f_before) return true;

        A_[idx(p, q)] = old_pq;
        A_[idx(q, p)] = old_qp;
        H_(r1, r2) = std::polar(1.0, TAU * old_pq);
        H_(r2, r1) = std::polar(1.0, TAU * old_qp);
        refresh_row(r1);
        if (r2 != r1) refresh_row(r2);
        F_ = f_before;
        return false;
    }

    void set_mirrored(int p, int q, double phase) {
        double a_pq = phase, a_qp = phase;
        if (opt_.symmetry == Symmetry::Hermitian) {
            if (p == q) { a_pq = a_qp = (phase < 0.5) ? 0.0 : 0.5; }
            else        { a_qp = wrap01(-phase); }
        }
        A_[idx(p, q)] = a_pq;
        A_[idx(q, p)] = a_qp;
        H_(p + 1, q + 1) = std::polar(1.0, TAU * a_pq);
        H_(q + 1, p + 1) = std::polar(1.0, TAU * a_qp);
    }

    // ---- exact single-coordinate minimiser ----------------------------------
    // With u = H_rc on the unit circle and C_rb the part of G_rb independent of
    // u, the r-part of F reads   const + 4 Re(u * Z),  Z = sum_{b != r} conj(C_rb) conj(H_bc).
    // Re(e^{i t} Z) is minimal at  u = -conj(Z)/|Z|.
    double exact_min_phase(int r, int c) const {
        const cd u = H_(r, c);
        cd Z(0.0, 0.0);
        for (int b = 0; b < n_; ++b) {
            if (b == r) continue;
            const cd w = std::conj(H_(b, c));   // |w| = 1
            const cd C = G_(r, b) - u * w;      // drop the k = c term
            Z += std::conj(C) * w;
        }
        const double az = std::abs(Z);
        if (!(az > 1e-300)) return phase01(u);
        return phase01(-std::conj(Z) / az);
    }

    const GenOptions& opt_;
    std::mt19937_64&  rng_;
    int n_, m_;
    std::vector<double> A_;
    CMat   H_, G_;
    double F_ = 0.0;
};

}  // namespace

bool hermitian_chm_possible(int n) {
    if (n <= 0) return false;
    if (n % 2 == 0) return true;
    const int r = int(std::lround(std::sqrt(double(n))));
    return r * r == n;
}

GenResult generate(const GenOptions& opt_in) {
    GenOptions opt = opt_in;
    GenResult  res;

    if (opt.size < 1) return res;
    if (opt.seed == 0) {
        std::random_device rd;
        opt.seed = (std::uint64_t(rd()) << 32) ^ std::uint64_t(rd()) ^
                   std::uint64_t(std::time(nullptr));
    }
    res.seed = opt.seed;

    if (opt.size == 1) {
        res.H = CMat::identity(1);
        res.H(0, 0) = cd(1.0, 0.0);
        res.ok = true;
        return res;
    }

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = (opt.threads > 0) ? opt.threads : omp_get_max_threads();
    if (nthreads < 1) nthreads = 1;
#endif
    // Independent restarts are embarrassingly parallel.  Pin BLAS to one thread
    // inside the region so the two thread pools do not oversubscribe the cores.
    blas_threads(1);

    std::atomic<bool> found{false};
    CMat              winner;
    long long         win_iters = 0;
    double            win_n1 = 0.0, win_nh = 0.0;
    long long         total_restarts = 0;
    double            global_best = std::numeric_limits<double>::infinity();
    double            best_n1 = 0.0, best_nh = 0.0;
    const auto        t_start = std::chrono::steady_clock::now();

    // Restart budget per thread (0 = unlimited: keep trying until success).
    const long long budget =
        (opt.restarts > 0) ? std::max<long long>(1, (opt.restarts + nthreads - 1) / nthreads) : 0;

#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads) reduction(+ : total_restarts)
#endif
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::mt19937_64 rng(opt.seed + 0x9E3779B97F4A7C15ULL * std::uint64_t(tid + 1));
        long long my_restarts = 0;
        double    my_best = std::numeric_limits<double>::infinity();
        double    my_best_n1 = 0.0, my_best_nh = 0.0;

        while (!found.load(std::memory_order_relaxed)) {
            if (budget > 0 && my_restarts >= budget) break;
            if (opt.max_seconds > 0.0 &&
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start)
                        .count() > opt.max_seconds)
                break;
            ++my_restarts;

            CMat      X;
            long long iters = 0;
            double    v1 = 0.0, vh = 0.0;
            bool      ok = false;

            if (opt.method == Method::Sinkhorn) {
                ok = sinkhorn_restart(X, opt, rng, iters, v1, vh);
                // Alternating projection converges only linearly and, past
                // N ~ 12, typically parks at n1 ~ 1e-5 after the phase kicks.
                // Taking the phases of that near-miss and running exact
                // coordinate descent on nh^2 closes the remaining gap.
                if (!ok && opt.polish && X.n == opt.size && v1 < 1e-2) {
                    RWCP  pol(opt, rng);
                    double p1 = 0.0, ph = 0.0;
                    pol.seed_from(X);
                    if (pol.sharpen(std::max<long long>(1 << 17, 512LL * opt.size * opt.size),
                                    p1, ph)) {
                        X = pol.matrix();
                        v1 = p1;
                        vh = ph;
                        ok = true;
                    }
                }
            } else {
                RWCP walker(opt, rng);
                ok = walker.restart(iters, v1, vh);
                X  = walker.matrix();
            }

            if (ok) {
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    if (!found.load(std::memory_order_relaxed)) {
                        found.store(true, std::memory_order_relaxed);
                        winner    = X;
                        win_iters = iters;
                        win_n1    = v1;
                        win_nh    = vh;
                    }
                }
                break;
            }
            if (v1 + vh < my_best) { my_best = v1 + vh; my_best_n1 = v1; my_best_nh = vh; }
            if (opt.verbose) { std::fputc('.', stderr); std::fflush(stderr); }
        }
        total_restarts += my_restarts;
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            if (my_best < global_best) {
                global_best = my_best;
                best_n1     = my_best_n1;
                best_nh     = my_best_nh;
            }
        }
    }

    blas_threads(0);  // 0 = restore the library default

    if (opt.verbose) std::fputc('\n', stderr);

    res.ok         = found.load();
    res.restarts   = total_restarts;
    res.iterations = win_iters;
    if (found) {
        res.H      = dephase(winner);
        // Dephasing with this convention preserves both H = H^T and H = H^dag,
        // so the symmetry requested on the command line survives into the file.
        res.n1_val = n1(res.H);
        res.nh_val = nh(res.H);
    } else {
        // Report the closest approach seen across all threads, not the last one.
        res.n1_val = best_n1;
        res.nh_val = best_nh;
        (void)win_n1; (void)win_nh;
    }
    return res;
}

// -----------------------------------------------------------------------------
// I/O
// -----------------------------------------------------------------------------
std::string timestamp() {
    const std::time_t t = std::time(nullptr);
    std::tm tmv{};
    localtime_r(&t, &tmv);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y%m%dT%H%M%S", &tmv);
    return std::string(buf);
}

std::string free_path(const std::string& path) {
    {
        std::ifstream probe(path);
        if (!probe) return path;
    }
    std::string stem = path, ext;
    const std::size_t slash = path.find_last_of('/');
    const std::size_t dot   = path.find_last_of('.');
    if (dot != std::string::npos && (slash == std::string::npos || dot > slash)) {
        stem = path.substr(0, dot);
        ext  = path.substr(dot);
    }
    for (int k = 1; k < 10000; ++k) {
        const std::string cand = stem + "_" + std::to_string(k) + ext;
        std::ifstream probe(cand);
        if (!probe) return cand;
    }
    return path;
}

namespace {

// Pull the first bracketed numeric array out of an Octave source file.
// Accepts:  function A = name ... A = [ ... ];  end
//           A = [ ... ];
//           a bare whitespace/newline separated table of numbers (.data files)
bool parse_numeric_array(const std::string& raw, std::vector<std::vector<double>>& rows,
                         std::string& err) {
    // 1. strip comments (% and # to end of line) outside of brackets
    std::string s;
    s.reserve(raw.size());
    bool in_comment = false;
    for (char c : raw) {
        if (in_comment) { if (c == '\n') { in_comment = false; s.push_back(c); } continue; }
        if (c == '%' || c == '#') { in_comment = true; continue; }
        s.push_back(c);
    }

    // 2. locate the array body
    std::string body;
    const std::size_t lb = s.find('[');
    if (lb != std::string::npos) {
        int depth = 0;
        std::size_t i = lb, rb = std::string::npos;
        for (; i < s.size(); ++i) {
            if (s[i] == '[') ++depth;
            else if (s[i] == ']') { if (--depth == 0) { rb = i; break; } }
        }
        if (rb == std::string::npos) { err = "unbalanced '[' in file"; return false; }
        body = s.substr(lb + 1, rb - lb - 1);
    } else {
        // no brackets: treat the whole (comment-stripped) file as a table, but
        // drop any line that carries non-numeric syntax such as "function"/"end"
        std::istringstream in(s);
        std::string line;
        while (std::getline(in, line)) {
            bool numeric_only = true;
            for (char c : line)
                if (!(std::isdigit(static_cast<unsigned char>(c)) || std::isspace(static_cast<unsigned char>(c)) ||
                      c == '+' || c == '-' || c == '.' || c == 'e' || c == 'E' || c == ',' || c == ';')) {
                    numeric_only = false; break;
                }
            if (numeric_only) { body += line; body += '\n'; }
        }
    }

    // 3. honour '...' continuations, then split rows on ';' and newlines
    for (std::size_t p = body.find("..."); p != std::string::npos; p = body.find("...", p)) {
        std::size_t e = p + 3;
        while (e < body.size() && body[e] != '\n') ++e;
        if (e < body.size()) ++e;
        body.replace(p, e - p, " ");
    }
    for (char& c : body) if (c == ',') c = ' ';

    rows.clear();
    std::string cur;
    auto flush_row = [&]() {
        std::istringstream ls(cur);
        std::vector<double> r;
        double v;
        while (ls >> v) r.push_back(v);
        if (!r.empty()) rows.push_back(std::move(r));
        cur.clear();
    };
    for (char c : body) {
        if (c == ';' || c == '\n') flush_row();
        else cur.push_back(c);
    }
    flush_row();

    if (rows.empty()) { err = "no numeric data found"; return false; }
    const std::size_t w = rows[0].size();
    for (const auto& r : rows)
        if (r.size() != w) { err = "ragged array: rows of unequal length"; return false; }
    if (w != rows.size()) {
        char b[128];
        std::snprintf(b, sizeof(b), "array is %zux%zu, expected a square array", rows.size(), w);
        err = b;
        return false;
    }
    return true;
}

}  // namespace

LoadResult load_phase_file(const std::string& path, Layout layout) {
    LoadResult r;
    std::ifstream f(path, std::ios::binary);
    if (!f) { r.error = "cannot open '" + path + "'"; return r; }
    std::string raw((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

    std::vector<std::vector<double>> A;
    if (!parse_numeric_array(raw, A, r.error)) return r;
    r.array_size = int(A.size());

    // Phases are documented to live in [0,1).  Accept radians too, so that
    // hand-made files are not silently misread: if anything exceeds 1 in
    // modulus but everything fits in [-2pi, 2pi], rescale by 1/(2pi).
    double amax = 0.0;
    for (const auto& row : A) for (double v : row) amax = std::max(amax, std::abs(v));
    if (amax > 1.0 + 1e-9 && amax <= TAU + 1e-9)
        for (auto& row : A) for (double& v : row) v /= TAU;

    const CMat core = from_core_phases(A);
    const CMat full = from_full_phases(A);

    if (layout == Layout::Core) { r.H = core; r.used = Layout::Core; r.ok = true; return r; }
    if (layout == Layout::Full) { r.H = full; r.used = Layout::Full; r.ok = true; return r; }

    const double ec = nh(core) + n1(core);
    const double ef = nh(full) + n1(full);
    if (ec <= ef) { r.H = core; r.used = Layout::Core; } else { r.H = full; r.used = Layout::Full; }
    r.ok = true;
    return r;
}

bool write_core_m(const std::string& path, const CMat& H, const std::string& stamp,
                  const std::string& provenance, int digits) {
    if (digits < 6)  digits = 6;
    if (digits > 17) digits = 17;
    // Octave requires the function name to match the file's base name, so the
    // identifier is derived from the path rather than fixed to "core_phases".
    std::string base = path;
    const std::size_t slash = base.find_last_of('/');
    if (slash != std::string::npos) base = base.substr(slash + 1);
    const std::size_t dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0, dot);
    std::string fname;
    for (std::size_t i = 0; i < base.size(); ++i)
        fname.push_back(is_ident_char(base[i], i == 0) ? base[i] : '_');
    if (fname.empty() || !is_ident_char(fname[0], true)) fname = "m_" + fname;

    const auto A = core_phases(H);

    std::ofstream f(path);
    if (!f) return false;
    f << "function A = " << fname << "\n";
    f << "% ---------------------------------------------------------------------"
         "---------\n";
    f << "% " << stamp << "\n";
    f << "% core phases\n";
    f << "% usage: >> A = " << fname << "; H = [ones(1, columns(A) + 1); "
      << "[ones(rows(A), 1), exp(2j*pi*A)]];\n";
    f << "% phases are normalised to the unit interval: A(j,k) in [0, 1)\n";
    if (!provenance.empty()) {
        std::istringstream ps(provenance);
        std::string line;
        while (std::getline(ps, line)) f << "% " << line << "\n";
    }
    f << "% ---------------------------------------------------------------------"
         "---------\n";
    f << "    A = [\n";
    for (std::size_t i = 0; i < A.size(); ++i) {
        f << "       ";
        for (std::size_t j = 0; j < A[i].size(); ++j) {
            char buf[64];
            std::snprintf(buf, sizeof(buf), " %.*g", digits, A[i][j]);
            f << buf;
        }
        f << (i + 1 < A.size() ? ";\n" : "\n");
    }
    f << "    ];\n";
    f << "end\n";
    return bool(f);
}

bool write_lambda_data(const std::string& path, const std::vector<cd>& lambda) {
    std::ofstream f(path);
    if (!f) return false;
    for (cd z : lambda) {
        char buf[96];
        std::snprintf(buf, sizeof(buf), "%.15f\t%.15f\n", z.real(), z.imag());
        f << buf;
    }
    return bool(f);
}

// =============================================================================
// Naming a result file after the matrix it holds
// =============================================================================

std::string exe_dir(const std::string& argv0) {
    char buf[4096];
    const ssize_t k = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (k > 0) {
        std::string p(buf, static_cast<std::size_t>(k));
        const std::size_t s = p.find_last_of('/');
        if (s != std::string::npos) return p.substr(0, s + 1);
    }
    const std::size_t slash = argv0.find_last_of('/');
    return (slash == std::string::npos) ? std::string("./")
                                        : argv0.substr(0, slash + 1);
}

std::string shell_quote(const std::string& s) {
    std::string q = "'";
    for (char c : s) {
        if (c == '\'') q += "'\\''";
        else           q += c;
    }
    return q + "'";
}

bool run_int(const std::string& cmd, long long& out) {
    std::FILE* p = popen((cmd + " 2>/dev/null").c_str(), "r");
    if (!p) return false;
    char buf[256] = {0};
    const bool got = (std::fgets(buf, sizeof(buf), p) != nullptr);
    const int status = pclose(p);
    if (!got || status != 0) return false;
    char* end = nullptr;
    const long long v = std::strtoll(buf, &end, 10);
    if (end == buf) return false;
    while (end && (*end == ' ' || *end == '\n' || *end == '\r' || *end == '\t')) ++end;
    if (end && *end != '\0') return false;
    out = v;
    return true;
}

char symmetry_letter(const CMat& H, double tol) {
    if (sym_error(H)  <= tol) return 'S';
    if (herm_error(H) <= tol) return 'H';
    return 'Y';
}

std::string name_after_matrix(const std::string& path, const CMat& H,
                              const std::string& stamp, const std::string& provenance,
                              int digits, const std::string& bindir,
                              std::FILE* log, long long d, long long L) {
    const char x = symmetry_letter(H);

    bool okd = (d >= 0), okL = (L >= 0);
    if (!okd) {
        if (log) { std::fprintf(log, "* now d (the defect) is being calculated ...\n");
                   std::fflush(log); }
        okd = run_int(shell_quote(bindir + "defect") + " --file " + shell_quote(path) +
                          " --bare", d);
        if (log) {
            if (okd) std::fprintf(log, "*   d = %lld\n", d);
            else     std::fprintf(log, "*   could not run %sdefect\n", bindir.c_str());
        }
    }
    if (!okL) {
        if (log) { std::fprintf(log, "* now L (the Haagerup invariants) is being calculated ...\n");
                   std::fflush(log); }
        okL = run_int(shell_quote(bindir + "lambda") + " --file " + shell_quote(path) +
                          " --bare --no-write", L);
        if (log) {
            if (okL) std::fprintf(log, "*   L = %lld\n", L);
            else     std::fprintf(log, "*   could not run %slambda\n", bindir.c_str());
        }
    }
    if (!okd || !okL) {
        if (log) std::fprintf(log, "*   keeping the neutral name %s\n", path.c_str());
        return path;
    }

    char nm[256];
    std::snprintf(nm, sizeof(nm), "%cH_%d_%lld_%lld_%s.m", x, H.n, d, L, stamp.c_str());
    const std::string dst = free_path(nm);

    // Not std::rename: Octave resolves a function by the file's base name, so
    // the "function A = ..." line has to follow the file.  Re-serialising from
    // the matrix in memory writes it correctly and loses nothing.
    if (!write_core_m(dst, H, stamp, provenance, digits)) {
        if (log) std::fprintf(log, "*   cannot write %s, keeping %s\n",
                              dst.c_str(), path.c_str());
        return path;
    }
    std::remove(path.c_str());
    return dst;
}

bool parse_matrix_name(const std::string& path, int& N, long long& d, long long& L) {
    std::string b = path;
    const std::size_t slash = b.find_last_of('/');
    if (slash != std::string::npos) b = b.substr(slash + 1);

    if (b.size() < 3 || b[1] != 'H' || b[2] != '_') return false;
    if (b[0] != 'S' && b[0] != 'H' && b[0] != 'Y')  return false;

    // three unsigned integers separated by '_', then the timestamp
    long long v[3];
    std::size_t i = 3;
    for (int f = 0; f < 3; ++f) {
        const std::size_t start = i;
        while (i < b.size() && std::isdigit(static_cast<unsigned char>(b[i]))) ++i;
        if (i == start || i >= b.size() || b[i] != '_') return false;
        v[f] = std::strtoll(b.substr(start, i - start).c_str(), nullptr, 10);
        ++i;                                  // step over the '_'
    }
    if (v[0] < 1 || v[0] > 100000) return false;
    N = static_cast<int>(v[0]);
    d = v[1];
    L = v[2];
    return true;
}

}  // namespace chm
