// =============================================================================
// check_me.cpp -- check MONOMIAL EQUIVALENCE of two complex Hadamard matrices
//
// 2026-09-04  W. Bruzda  <w.bruzda@cft.edu.pl>
//
// Two CHM are monomially equivalent when
//
//     H2 = D1 P1 H1 P2 D2
//
// with D1, D2 diagonal unitary and P1, P2 permutations.  Dephasing H2 at one of
// its entries (a,b),
//
//     dephase_{a,b}(H)_{kl} = H_kl H_ab / (H_al H_kb) ,
//
// turns row a and column b into ones and so fixes the two diagonal factors
// completely; only the permutations are left.  Every dephased representative of
// the class of H2 is therefore
//
//     P_row . dephase_{a,b}(H2) . P_col
//
// for some base (a,b).  The program walks all N*N bases and, for each, searches
// for the permutations explicitly: rows are assigned one at a time and a branch
// is cut as soon as the columns of the partial matrix can no longer be matched.
// The search is exhaustive, so a negative answer is a proof (up to the phase
// tolerance), and a positive answer comes with the witnesses sigma and tau,
// which are then verified on the actual phases.
//
// C++ port of find_equiv.py, reading the .m files written by get_chm directly:
//
//     ./check_me SH_0_139_20260904T162711.m SH_0_139_20260905T081530.m
//
// Build:
//   g++ -std=c++17 -O3 -march=znver3 -mtune=znver3 -fopenmp -flto
//       -o check_me check_me.cpp chm.cpp -llapacke -llapack -lopenblas -lm
// =============================================================================

#define CHM_PROGRAM_ABI 20260908   // must match CHM_ABI in chm.hpp
#include "chm.hpp"
#if CHM_ABI < CHM_PROGRAM_ABI
#error "MIXED SOURCE TREE: chm.hpp is older than this .cpp. Delete the directory and extract one archive into it; see COMPILE.txt."
#endif

#include "cli.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace chm;

namespace {

using Phase = std::vector<std::vector<double>>;   // N x N phases in [0,1)

inline double wrap01(double x) {
    x = std::fmod(x, 1.0);
    if (x < 0.0) x += 1.0;
    return x;
}

// N x N phase matrix of H (entries are unimodular, so only the phase matters)
Phase phases_of(const CMat& H) {
    const int n = H.n;
    Phase A(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) A[i][j] = wrap01(std::arg(H(i, j)) / TAU);
    return A;
}

Phase dephase_at(const Phase& A, int a, int b) {
    const int n = int(A.size());
    Phase D(n, std::vector<double>(n, 0.0));
    for (int k = 0; k < n; ++k)
        for (int l = 0; l < n; ++l)
            D[k][l] = wrap01(A[k][l] + A[a][b] - A[a][l] - A[k][b]);
    return D;
}

// the (N-1)x(N-1) core: drop row a and column b, which are now all ones
Phase core_of(const Phase& D, int a, int b) {
    const int n = int(D.size());
    Phase C;
    for (int k = 0; k < n; ++k) {
        if (k == a) continue;
        std::vector<double> row;
        for (int l = 0; l < n; ++l)
            if (l != b) row.push_back(D[k][l]);
        C.push_back(std::move(row));
    }
    return C;
}

// ---------------------------------------------------------------------------
// Codebook: map phases to small integers so the permutation search is exact.
// Built from every phase that any dephasing of either matrix can produce.
// ---------------------------------------------------------------------------
class Codebook {
  public:
    explicit Codebook(double tol) : tol_(tol) {}

    void add(double v) { raw_.push_back(wrap01(v)); }

    void build() {
        std::sort(raw_.begin(), raw_.end());
        for (double v : raw_)
            if (reps_.empty() || v - reps_.back() >= tol_) reps_.push_back(v);
        // 0.999... and 0.000... are the same phase
        if (reps_.size() > 1 && (1.0 - reps_.back()) + reps_.front() < tol_)
            reps_.pop_back();
        raw_.clear();
        raw_.shrink_to_fit();
    }

    int label(double v) const {
        v = wrap01(v);
        auto it = std::lower_bound(reps_.begin(), reps_.end(), v);
        int best = -1;
        double bd = 2.0;
        const int m = int(reps_.size());
        const int p = int(it - reps_.begin());
        for (int c : {p - 1, p, p + 1, 0, m - 1}) {
            if (c < 0 || c >= m) continue;
            double d = std::fabs(reps_[std::size_t(c)] - v);
            d = std::min(d, 1.0 - d);                    // circular distance
            if (d < bd) { bd = d; best = c; }
        }
        return best;
    }

    std::size_t size() const { return reps_.size(); }

    double min_gap() const {
        if (reps_.size() < 2) return 1.0;
        double g = 1.0;
        for (std::size_t i = 1; i < reps_.size(); ++i)
            g = std::min(g, reps_[i] - reps_[i - 1]);
        return std::min(g, (1.0 - reps_.back()) + reps_.front());
    }

  private:
    double tol_;
    std::vector<double> raw_;
    std::vector<double> reps_;
};

using IMat = std::vector<std::vector<int>>;

IMat labelled(const Phase& C, const Codebook& book) {
    IMat L(C.size(), std::vector<int>(C.size(), 0));
    for (std::size_t i = 0; i < C.size(); ++i)
        for (std::size_t j = 0; j < C.size(); ++j) L[i][j] = book.label(C[i][j]);
    return L;
}

// ---------------------------------------------------------------------------
// Row assignment with column-prefix pruning.
// After choosing sigma(0..depth-1), the multiset of column prefixes of the two
// partial matrices must agree; otherwise no completion can exist.
// ---------------------------------------------------------------------------
struct RowSearch {
    const IMat& D;
    const IMat& T;
    int m;
    long long nodes = 0;
    std::vector<int> sigma;
    std::vector<char> used;

    RowSearch(const IMat& d, const IMat& t, int mm)
        : D(d), T(t), m(mm), used(std::size_t(mm), 0) {}

    bool feasible() const {
        const int depth = int(sigma.size());
        std::vector<std::vector<int>> ct(static_cast<std::size_t>(m));
        std::vector<std::vector<int>> cd(static_cast<std::size_t>(m));
        for (int j = 0; j < m; ++j) {
            ct[std::size_t(j)].reserve(std::size_t(depth));
            cd[std::size_t(j)].reserve(std::size_t(depth));
            for (int i = 0; i < depth; ++i) {
                ct[std::size_t(j)].push_back(T[std::size_t(i)][std::size_t(j)]);
                cd[std::size_t(j)].push_back(D[std::size_t(sigma[std::size_t(i)])][std::size_t(j)]);
            }
        }
        std::sort(ct.begin(), ct.end());
        std::sort(cd.begin(), cd.end());
        return ct == cd;
    }

    bool rec() {
        ++nodes;
        if (int(sigma.size()) == m) return true;
        for (int r = 0; r < m; ++r) {
            if (used[std::size_t(r)]) continue;
            used[std::size_t(r)] = 1;
            sigma.push_back(r);
            if (feasible() && rec()) return true;
            sigma.pop_back();
            used[std::size_t(r)] = 0;
        }
        return false;
    }

    // Same walk, but collecting every complete assignment instead of stopping
    // at the first.  Used by --count, which needs all of them.
    std::vector<std::vector<int>> all;
    void rec_all() {
        ++nodes;
        if (int(sigma.size()) == m) { all.push_back(sigma); return; }
        for (int r = 0; r < m; ++r) {
            if (used[std::size_t(r)]) continue;
            used[std::size_t(r)] = 1;
            sigma.push_back(r);
            if (feasible()) rec_all();
            sigma.pop_back();
            used[std::size_t(r)] = 0;
        }
    }
};

// Once the rows are fixed, matching the columns is exact vector equality.
bool column_perm(const IMat& D, const IMat& T, const std::vector<int>& sigma,
                 int m, std::vector<int>& tau) {
    tau.clear();
    std::vector<char> used(std::size_t(m), 0);
    for (int j = 0; j < m; ++j) {
        int hit = -1;
        for (int c = 0; c < m && hit < 0; ++c) {
            if (used[std::size_t(c)]) continue;
            bool ok = true;
            for (int i = 0; i < m && ok; ++i)
                ok = D[std::size_t(sigma[std::size_t(i)])][std::size_t(c)] ==
                     T[std::size_t(i)][std::size_t(j)];
            if (ok) hit = c;
        }
        if (hit < 0) return false;
        used[std::size_t(hit)] = 1;
        tau.push_back(hit);
    }
    return true;
}

void usage(const char* prog) {
    std::printf(
        "usage: %s FILE1 FILE2 [--tol T] [--quiet]\n"
        "\n"
        "  FILE1 FILE2   core-phase files as written by get_chm, e.g.\n"
        "                SH_0_139_20260904T162711.m   (.hp files also work)\n"
        "  --tol T       two phases count as equal below T      (default 1e-9)\n"
        "                the files hold ~17 digits, so 1e-9 is right for them;\n"
        "                lower it only for refined high-precision input\n"
        "  --count       count every monomial equivalence FILE2 -> FILE1 instead\n"
        "                of stopping at the first.  With the same file twice this\n"
        "                is the order of the automorphism group modulo scalars,\n"
        "                |Aut(H)/U(1)| -- the group that governs whether a\n"
        "                self-transpose class has a symmetric representative\n"
        "  --quiet       print the verdict (or, with --count, the number) alone\n"
        "\n"
        "exit: 0 = monomially equivalent, 1 = not equivalent, 2 = bad usage\n",
        prog);
}

}  // namespace

int main(int argc, char** argv) {
    cli::Args args = cli::parse(argc, argv);

    if (args.has("--help") || args.has("-h")) { usage(args.prog.c_str()); return 0; }

    // Split the command line into flags and the two positional file names.
    // Only --tol takes a value, so the token after it is skipped.
    std::vector<std::string> files;
    for (std::size_t i = 0; i < args.tok.size(); ++i) {
        const std::string& t = args.tok[i];
        if (t == "--tol") { ++i; continue; }
        if (t == "--quiet" || t == "--count" || t == "--help" || t == "-h") continue;
        if (t.rfind("-", 0) == 0) {
            std::fprintf(stderr, "%s: unknown flag '%s'\n", args.prog.c_str(), t.c_str());
            return 2;
        }
        files.push_back(t);
    }
    if (files.size() != 2) {
        std::fprintf(stderr, "* warning: exactly two matrix files are required\n");
        usage(args.prog.c_str());
        return 2;
    }

    double tol = 1e-9;
    args.real("--tol", tol);
    if (!args.error.empty()) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), args.error.c_str());
        return 2;
    }
    const bool quiet = args.has("--quiet");
    const bool count = args.has("--count");

    LoadResult L1 = load_phase_file(files[0], Layout::Core);
    LoadResult L2 = load_phase_file(files[1], Layout::Core);
    if (!L1.ok) { std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), L1.error.c_str()); return 2; }
    if (!L2.ok) { std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), L2.error.c_str()); return 2; }

    if (L1.H.n != L2.H.n) {
        if (!quiet)
            std::printf("* different orders (%d and %d)\n", L1.H.n, L2.H.n);
        std::printf("* INEQUIVALENT\n");
        return 1;
    }
    const int n = L1.H.n;
    const int m = n - 1;

    const Phase A1 = phases_of(L1.H);
    const Phase A2 = phases_of(L2.H);

    if (!quiet) {
        std::printf("* %s   N = %d, nh = %.3e, n1 = %.3e\n", files[0].c_str(), n,
                    nh(L1.H), n1(L1.H));
        std::printf("* %s   N = %d, nh = %.3e, n1 = %.3e\n", files[1].c_str(), n,
                    nh(L2.H), n1(L2.H));
    }

    // one codebook over every phase any dephasing can produce, so that all the
    // comparisons below are exact integer comparisons
    Codebook book(tol);
    for (const Phase* A : {&A1, &A2})
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b) {
                const Phase D = dephase_at(*A, a, b);
                for (int k = 0; k < n; ++k)
                    for (int l = 0; l < n; ++l) book.add(D[k][l]);
            }
    book.build();
    if (!quiet) {
        std::printf("* codebook: %zu distinct phases, smallest gap %.3e (tolerance %g)\n",
                    book.size(), book.min_gap(), tol);
        if (book.min_gap() < 100 * tol)
            std::printf("  WARNING: phases lie closer together than 100x the tolerance;\n"
                        "           refine the input or lower --tol before trusting this\n");
    }

    const IMat T = labelled(core_of(dephase_at(A1, 0, 0), 0, 0), book);

    // ---- --count: enumerate every equivalence instead of the first ----------
    if (count) {
        long long total = 0, nodes_c = 0;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b) {
                const IMat D = labelled(core_of(dephase_at(A2, a, b), a, b), book);
                RowSearch rs(D, T, m);
                rs.rec_all();
                nodes_c += rs.nodes;
                // the columns of a CHM are pairwise non-proportional, so for a
                // given base and row assignment the column permutation, if it
                // exists at all, is unique -- hence one equivalence per hit
                for (const auto& sg : rs.all) {
                    std::vector<int> tau;
                    if (column_perm(D, T, sg, m, tau)) ++total;
                }
            }
        if (quiet) std::printf("%lld\n", total);
        else {
            const bool self = (files[0] == files[1]);
            std::printf("* %lld monomial %s (mod scalars), %lld search nodes\n",
                        total, self ? "automorphisms" : "equivalences", nodes_c);
            if (self)
                std::printf("* |Aut(H)/U(1)| = %lld%s\n", total,
                            total == 1 ? "   (trivial: only the scalars)" : "");
            std::printf(total > 0 ? "* EQUIVALENT\n" : "* INEQUIVALENT\n");
        }
        return total > 0 ? 0 : 1;
    }

    long long nodes = 0;
    for (int a = 0; a < n; ++a) {
        for (int b = 0; b < n; ++b) {
            const IMat D = labelled(core_of(dephase_at(A2, a, b), a, b), book);
            RowSearch rs(D, T, m);
            const bool got = rs.rec();
            nodes += rs.nodes;
            if (!got) continue;
            std::vector<int> tau;
            if (!column_perm(D, T, rs.sigma, m, tau)) continue;

            // verify on the actual phases, not on the labels
            const Phase Dp = core_of(dephase_at(A2, a, b), a, b);
            const Phase Tp = core_of(dephase_at(A1, 0, 0), 0, 0);
            double err = 0.0;
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < m; ++j) {
                    double d = std::fabs(Dp[std::size_t(rs.sigma[std::size_t(i)])]
                                            [std::size_t(tau[std::size_t(j)])] -
                                         Tp[std::size_t(i)][std::size_t(j)]);
                    err = std::max(err, std::min(d, 1.0 - d));
                }
            if (!quiet) {
                std::printf("* base of %s: (row %d, col %d)\n", files[1].c_str(), a, b);
                std::printf("* row permutation:   ");
                for (int v : rs.sigma) std::printf(" %d", v);
                std::printf("\n* column permutation:");
                for (int v : tau) std::printf(" %d", v);
                std::printf("\n* max phase mismatch after applying them: %.3e\n", err);
                std::printf("* (%lld search nodes)\n", nodes);
            }
            std::printf("* EQUIVALENT\n");
            return 0;
        }
    }

    if (!quiet)
        std::printf("* all %d bases exhausted, %lld search nodes\n", n * n, nodes);
    std::printf("* INEQUIVALENT\n");
    return 1;
}
