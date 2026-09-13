// =============================================================================
// symmetrize.cpp -- put a complex Hadamard matrix into symmetric form, if its
//                   monomial equivalence class admits one
//
// 2026-09-05  W. Bruzda  <w.bruzda@cft.edu.pl>
//
//     ./symmetrize --file xH_N_d_L_YYYYMMDDThhmmss.m
//
// -----------------------------------------------------------------------------
// WHAT IS BEING DECIDED
// -----------------------------------------------------------------------------
// Two CHM are monomially equivalent when  H2 = D1 P1 H1 P2 D2  with D1, D2
// diagonal unitary and P1, P2 permutations.  The question here is whether the
// class of the given H contains a representative S with S = S^T, and, if so, to
// produce one.  The search below is exhaustive, so "no" is a proof (up to the
// phase tolerance) and not a failure to find.
//
// -----------------------------------------------------------------------------
// THE REDUCTION THAT MAKES IT CHEAP
// -----------------------------------------------------------------------------
// Dephasing at an entry (a,b),
//
//     K = dephase_{a,b}(H),      K_kl = H_kl H_ab / (H_al H_kb) ,
//
// turns row a and column b into ones and thereby CANCELS both diagonal factors:
// for H' = D1 P H Q D2 one checks directly that
//
//     dephase_{a,b}(H') = P . dephase_{p(a), q(b)}(H) . Q                  (*)
//
// -- the d's drop out.  So every dephased member of the class of H is a row and
// column permutation of one of the N*N matrices dephase_{a,b}(H).
//
// Now suppose S = S^T lies in the class.  Dephasing S at a DIAGONAL entry (a,a)
// keeps it symmetric, because S_al = S_la and S_ka = S_ak give
//
//     dephase_{a,a}(S)_kl = S_kl S_aa / (S_al S_ka) = dephase_{a,a}(S)_lk .
//
// By (*) that symmetric matrix is P K Q for some base and some P, Q.  Write
// P K Q as (i,j) -> K_{alpha(i), beta(j)} and put rho = alpha^{-1} beta.  Then
//
//     symmetry  <=>  K_{alpha(i), alpha(rho(j))} = K_{alpha(j), alpha(rho(i))}
//
// and the matrix C_ij = K_{alpha(i), alpha(rho(j))} is a SIMULTANEOUS row and
// column relabelling of K' with K'_uv = K_{u, (alpha rho alpha^{-1})(v)}.  Such
// a relabelling preserves symmetry, and alpha rho alpha^{-1} runs over all
// permutations as rho does.  Hence
//
//     the class of H contains a symmetric matrix
//         <=>  for some base (a,b) and some permutation pi,
//              K_{i, pi(j)} = K_{j, pi(i)}  for all i, j,        K = dephase_{a,b}(H)
//
// i.e. ONE permutation has to be found per base, not two.  The witness is then
// S = K P_pi, already symmetric and already a CHM.
//
// Two further facts prune the search hard:
//
//   * pi(a) = b is forced.  Row a of K P_pi is all ones (row a of K is), so by
//     symmetry column a must be all ones too; column a of K P_pi is column
//     pi(a) of K, and column b is the only all-ones column of K -- a second one
//     would make two columns of H proportional, which orthogonality forbids.
//     So the symmetric representative is automatically dephased at (a,a).
//
//   * the multiset of row i of K must equal the multiset of column pi(i) of K,
//     since row i of K P_pi is a permutation of row i of K and has to equal
//     column i of K P_pi, which is column pi(i) of K.  Matching multisets is
//     checked once per base and usually leaves very few candidates per row.
//
// The remaining assignment pi(0), pi(1), ... is a backtracking search: fixing
// pi(j) = c requires K_{i,c} = K_{j,pi(i)} for every i already assigned, which
// is O(assigned) per candidate and cuts almost immediately.
//
// -----------------------------------------------------------------------------
// A NEGATIVE ANSWER
// -----------------------------------------------------------------------------
// "not symmetrisable" is a statement about THIS matrix, not about the family it
// may sit in: a matrix of positive defect belongs to a continuous family whose
// other members are inequivalent to it, and one of those may well be
// symmetrisable when this one is not.  For defect 0 the point is isolated and
// the statement is about the class outright.
//
// Build:
//   g++ -std=c++17 -O3 -march=znver3 -mtune=znver3 -fopenmp -flto
//       -o symmetrize symmetrize.cpp chm.cpp -llapacke -llapack -lopenblas -lm
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

// N x N phase matrix of H (every entry is unimodular, so only the phase matters)
Phase phases_of(const CMat& H) {
    const int n = H.n;
    Phase P(static_cast<std::size_t>(n), std::vector<double>(static_cast<std::size_t>(n)));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            P[std::size_t(i)][std::size_t(j)] =
                wrap01(std::arg(H(i, j)) / TAU);
    return P;
}

// dephase_{a,b}(H) in phases:  P_kl + P_ab - P_al - P_kb   (mod 1)
Phase dephase_at(const Phase& P, int a, int b) {
    const std::size_t n = P.size();
    Phase K(n, std::vector<double>(n));
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t l = 0; l < n; ++l)
            K[k][l] = wrap01(P[k][l] + P[std::size_t(a)][std::size_t(b)]
                             - P[std::size_t(a)][l] - P[k][std::size_t(b)]);
    return K;
}

// -----------------------------------------------------------------------------
// Codebook: the distinct phases occurring anywhere, so that comparisons become
// integer comparisons and the tolerance is applied exactly once.
// -----------------------------------------------------------------------------
class Codebook {
public:
    void add(double v)               { vals_.push_back(wrap01(v)); }
    void build(double tol) {
        std::sort(vals_.begin(), vals_.end());
        reps_.clear();
        for (double v : vals_)
            if (reps_.empty() || v - reps_.back() >= tol) reps_.push_back(v);
        // 0 and 1 are the same phase
        if (reps_.size() > 1 && (1.0 - reps_.back()) + reps_.front() < tol)
            reps_.pop_back();
        gap_ = 1.0;
        for (std::size_t i = 1; i < reps_.size(); ++i)
            gap_ = std::min(gap_, reps_[i] - reps_[i - 1]);
        if (reps_.size() > 1)
            gap_ = std::min(gap_, (1.0 - reps_.back()) + reps_.front());
        vals_.clear();
        vals_.shrink_to_fit();
    }
    // index of the representative nearest to v, cyclically
    int label(double v) const {
        v = wrap01(v);
        auto it = std::lower_bound(reps_.begin(), reps_.end(), v);
        std::size_t hi = static_cast<std::size_t>(it - reps_.begin());
        std::size_t lo = (hi == 0) ? reps_.size() - 1 : hi - 1;
        if (hi == reps_.size()) hi = 0;
        const double dhi = std::min(std::fabs(reps_[hi] - v), 1.0 - std::fabs(reps_[hi] - v));
        const double dlo = std::min(std::fabs(reps_[lo] - v), 1.0 - std::fabs(reps_[lo] - v));
        return static_cast<int>(dhi <= dlo ? hi : lo);
    }
    std::size_t size()    const { return reps_.size(); }
    double      min_gap() const { return gap_; }

private:
    std::vector<double> vals_, reps_;
    double gap_ = 1.0;
};

using IMat = std::vector<std::vector<int>>;

IMat labelled(const Phase& K, const Codebook& cb) {
    const std::size_t n = K.size();
    IMat M(n, std::vector<int>(n));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            M[i][j] = cb.label(K[i][j]);
    return M;
}

// -----------------------------------------------------------------------------
// The per-base search:  find pi with  M[i][pi(j)] == M[j][pi(i)]  for all i, j.
// -----------------------------------------------------------------------------
struct SymSearch {
    const IMat& M;
    int n;
    int a, b;                       // the base; pi(a) = b is forced
    std::vector<int> pi, used;
    std::vector<std::vector<int>> cand;   // columns whose multiset matches row i
    long long nodes = 0;

    SymSearch(const IMat& M_, int a_, int b_)
        : M(M_), n(static_cast<int>(M_.size())), a(a_), b(b_),
          pi(static_cast<std::size_t>(n), -1),
          used(static_cast<std::size_t>(n), 0) {}

    // multiset of row i must equal multiset of column pi(i)
    bool build_candidates() {
        std::vector<std::vector<int>> rows(static_cast<std::size_t>(n)),
                                      cols(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            rows[std::size_t(i)].resize(static_cast<std::size_t>(n));
            cols[std::size_t(i)].resize(static_cast<std::size_t>(n));
            for (int j = 0; j < n; ++j) {
                rows[std::size_t(i)][std::size_t(j)] = M[std::size_t(i)][std::size_t(j)];
                cols[std::size_t(i)][std::size_t(j)] = M[std::size_t(j)][std::size_t(i)];
            }
            std::sort(rows[std::size_t(i)].begin(), rows[std::size_t(i)].end());
            std::sort(cols[std::size_t(i)].begin(), cols[std::size_t(i)].end());
        }
        cand.assign(static_cast<std::size_t>(n), {});
        for (int i = 0; i < n; ++i) {
            for (int c = 0; c < n; ++c)
                if (rows[std::size_t(i)] == cols[std::size_t(c)])
                    cand[std::size_t(i)].push_back(c);
            if (cand[std::size_t(i)].empty()) return false;
        }
        // pi(a) = b is forced by the all-ones row/column argument
        if (std::find(cand[std::size_t(a)].begin(), cand[std::size_t(a)].end(), b)
            == cand[std::size_t(a)].end())
            return false;
        cand[std::size_t(a)] = {b};
        return true;
    }

    // consistency of pi(j) = c against everything already assigned
    bool ok_with_assigned(int j, int c) const {
        for (int i = 0; i < n; ++i) {
            if (pi[std::size_t(i)] < 0) continue;
            if (M[std::size_t(i)][std::size_t(c)] !=
                M[std::size_t(j)][std::size_t(pi[std::size_t(i)])]) return false;
        }
        return true;
    }

    // assign in an order that starts from the forced value and then takes the
    // most constrained row first
    bool rec(const std::vector<int>& order, std::size_t k) {
        if (k == order.size()) return true;
        ++nodes;
        const int j = order[k];
        for (int c : cand[std::size_t(j)]) {
            if (used[std::size_t(c)]) continue;
            if (!ok_with_assigned(j, c)) continue;
            pi[std::size_t(j)] = c;
            used[std::size_t(c)] = 1;
            if (rec(order, k + 1)) return true;
            used[std::size_t(c)] = 0;
            pi[std::size_t(j)] = -1;
        }
        return false;
    }

    bool cut_by_multisets = false;      // rejected before any backtracking

    bool run() {
        if (!build_candidates()) { cut_by_multisets = true; return false; }
        std::vector<int> order(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) order[std::size_t(i)] = i;
        std::stable_sort(order.begin(), order.end(), [&](int p, int q) {
            return cand[std::size_t(p)].size() < cand[std::size_t(q)].size();
        });
        return rec(order, 0);
    }
};

void usage(const char* prog) {
    std::printf(
        "usage: %s --file FILE [--layout core|full|auto] [--tol T]\n"
        "                 [--out FILE] [--no-write] [--all] [--quiet] [--bare]\n"
        "\n"
        "Decides whether the monomial equivalence class of the given complex\n"
        "Hadamard matrix contains a symmetric representative S = S^T, and writes\n"
        "one out when it does.  The search is exhaustive: a negative answer is a\n"
        "proof up to the phase tolerance.\n"
        "\n"
        "  --file FILE     the matrix, as written by get_chm; mandatory\n"
        "  --layout L      core | full | auto                     (default core)\n"
        "  --tol T         two phases count as equal below T      (default 1e-9)\n"
        "  --out FILE      output path; the default is built from the matrix,\n"
        "                  SH_<N>_<d>_<L>_<stamp>.m, reusing d and L from the\n"
        "                  input file name when it carries them (both are class\n"
        "                  invariants) and running ./defect and ./lambda if not\n"
        "  --no-write      decide only, write nothing\n"
        "  --all           try every base instead of stopping at the first hit,\n"
        "                  and report how many of the N^2 bases work\n"
        "  --quiet         print the output file name only, nothing on failure\n"
        "  --bare          print 1 (symmetrisable) or 0, and nothing else\n"
        "  --help, -h      this message\n"
        "\n"
        "exit: 0 = a symmetric representative was found, 1 = the class has none,\n"
        "      2 = bad usage or the file cannot be read\n",
        prog);
}

}  // namespace

int main(int argc, char** argv) {
    cli::Args args = cli::parse(argc, argv);

    const std::vector<std::string> known = {
        "--file", "--layout", "--tol", "--out", "--no-write", "--all",
        "--quiet", "--bare", "--help", "-h"};
    const std::vector<std::string> valued = {"--file", "--layout", "--tol", "--out"};

    if (args.has("--help") || args.has("-h")) { usage(args.prog.c_str()); return 0; }
    if (!args.check_known(known, valued)) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), args.error.c_str());
        usage(args.prog.c_str());
        return 2;
    }

    std::string file;
    if (!args.value("--file", file) || file.empty()) {
        std::fprintf(stderr, "* warning: --file is mandatory\n");
        usage(args.prog.c_str());
        return 2;
    }

    Layout layout = Layout::Core;
    std::string lv;
    if (args.value("--layout", lv)) {
        if      (lv == "core") layout = Layout::Core;
        else if (lv == "full") layout = Layout::Full;
        else if (lv == "auto") layout = Layout::Auto;
        else { std::fprintf(stderr, "%s: --layout must be core, full or auto\n",
                            args.prog.c_str()); return 2; }
    }

    double tol = 1e-9;
    args.real("--tol", tol);
    if (!(tol > 0.0)) {
        std::fprintf(stderr, "%s: --tol must be positive\n", args.prog.c_str());
        return 2;
    }

    const bool quiet    = args.has("--quiet");
    const bool bare     = args.has("--bare");
    const bool no_write = args.has("--no-write");
    const bool all      = args.has("--all");
    const bool verbose  = !quiet && !bare;

    std::string out;
    args.value("--out", out);

    LoadResult in = load_phase_file(file, layout);
    if (!in.ok) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), in.error.c_str());
        return 2;
    }
    const CMat& H = in.H;
    const int n = H.n;

    const double nh0 = nh(H), n10 = n1(H);
    if (verbose)
        std::printf("* %s   N = %d, nh = %.3e, n1 = %.3e\n",
                    file.c_str(), n, nh0, n10);

    // The argument behind a negative answer uses orthogonality twice: to know
    // that column b is the only all-ones column of the dephased matrix (which
    // is what forces pi(a) = b), and to know that the witness really is a CHM.
    // On a matrix that is not one, both the verdict and the witness would be
    // meaningless, so this is refused rather than answered.
    if (nh0 > 1e-7 || n10 > 1e-7) {
        std::fflush(stdout);            // keep the two streams in order
        std::fprintf(stderr,
                     "%s: '%s' is not a complex Hadamard matrix "
                     "(nh = %.3e, n1 = %.3e).\n"
                     "  Symmetrisability is only defined on the monomial class of a CHM,\n"
                     "  and the search below assumes orthogonality; refusing to answer.\n"
                     "  Check --layout, or run ./summary --file '%s' to see what was read.\n",
                     args.prog.c_str(), file.c_str(), nh0, n10, file.c_str());
        return 2;
    }

    // ---- all N*N dephasings, one shared codebook -----------------------------
    const Phase P = phases_of(H);
    std::vector<Phase> K(static_cast<std::size_t>(n) * n);
    Codebook cb;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
            Phase Kab = dephase_at(P, a, b);
            for (const auto& row : Kab)
                for (double v : row) cb.add(v);
            K[std::size_t(a) * n + b] = std::move(Kab);
        }
    cb.build(tol);

    if (verbose)
        std::printf("* codebook: %zu distinct phases, smallest gap %.3e (tolerance %g)\n",
                    cb.size(), cb.min_gap(), tol);
    if (cb.min_gap() < 100.0 * tol && verbose)
        std::printf("*   warning: the smallest gap is close to the tolerance;\n"
                    "*   the phases may not be resolved -- refine the matrix or raise --tol\n");

    // ---- the search ---------------------------------------------------------
    bool found = false;
    int best_a = -1, best_b = -1;
    std::vector<int> best_pi;
    long long nodes = 0, hits = 0, cut = 0;

    for (int a = 0; a < n && (all || !found); ++a)
        for (int b = 0; b < n && (all || !found); ++b) {
            const IMat M = labelled(K[std::size_t(a) * n + b], cb);
            SymSearch s(M, a, b);
            const bool got = s.run();
            nodes += s.nodes;
            if (s.cut_by_multisets) ++cut;
            if (!got) continue;
            ++hits;
            if (!found) { found = true; best_a = a; best_b = b; best_pi = s.pi; }
        }

    if (bare) { std::printf("%d\n", found ? 1 : 0); return found ? 0 : 1; }

    if (!found) {
        if (!quiet)
            std::printf("* all %d bases exhausted: %lld ruled out by the row/column\n"
                        "*   multisets alone, the rest by %lld backtracking nodes\n"
                        "* NOT SYMMETRISABLE"
                        "  (no monomial transformation of this matrix is symmetric)\n",
                        n * n, cut, nodes);
        return 1;
    }

    // ---- build the witness  S = K P_pi  and check it ------------------------
    const Phase& Kb = K[std::size_t(best_a) * n + best_b];
    CMat S(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const double ph = Kb[std::size_t(i)][std::size_t(best_pi[std::size_t(j)])];
            S(i, j) = std::polar(1.0, TAU * ph);
        }

    const double se = sym_error(S), nhv = nh(S), n1v = n1(S);

    if (!quiet) {
        std::printf("* base: dephase at (row %d, col %d)\n", best_a, best_b);
        std::printf("* column permutation: ");
        for (int j = 0; j < n; ++j) std::printf("%d ", best_pi[std::size_t(j)]);
        std::printf("\n");
        std::printf("* || S - S^T ||_F = %.3e,  nh = %.3e,  n1 = %.3e\n", se, nhv, n1v);
        if (all)
            std::printf("* %lld of the %d bases admit a symmetric form (%lld nodes)\n",
                        hits, n * n, nodes);
        else
            std::printf("* (%lld search nodes)\n", nodes);
    }

    if (se > 1e-7 || nhv > 1e-7) {
        std::fprintf(stderr, "%s: the witness failed its own check "
                             "(|S-S^T| = %.3e, nh = %.3e); refusing to write it\n",
                     args.prog.c_str(), se, nhv);
        return 1;
    }

    if (no_write) {
        if (!quiet) std::printf("* SYMMETRISABLE   (--no-write: nothing written)\n");
        return 0;
    }

    // ---- write the symmetric representative ---------------------------------
    const std::string stamp = timestamp();
    char prov[512];
    std::snprintf(prov, sizeof(prov),
                  "N = %d, symmetric representative of the monomial class of\n"
                  "%s\n"
                  "obtained by dephasing at (%d, %d) and permuting columns\n"
                  "|| S - S^T ||_F = %.3e\n"
                  "nh = || H*H' - N*I ||_F = %.3e\n"
                  "n1 = ||  |H| - 1     ||_F = %.3e",
                  n, file.c_str(), best_a, best_b, se, nhv, n1v);

    const bool auto_name = out.empty();
    if (auto_name) out = free_path("A_" + stamp + ".m");

    if (!write_core_m(out, S, stamp, prov, 17)) {
        std::fprintf(stderr, "%s: cannot write '%s'\n", args.prog.c_str(), out.c_str());
        return 1;
    }

    if (auto_name) {
        // d and #Lambda are invariants of the class, so when the input file name
        // carries them there is nothing to recompute.
        int nn = 0; long long d = -1, L = -1;
        if (parse_matrix_name(file, nn, d, L) && nn == n) {
            if (!quiet)
                std::printf("* d = %lld and L = %lld taken from the input file name "
                            "(both are class invariants)\n", d, L);
        } else {
            d = L = -1;
        }
        out = name_after_matrix(out, S, stamp, prov, 17, exe_dir(args.prog),
                                quiet ? stderr : stdout, d, L);
    }

    if (quiet) std::printf("%s\n", out.c_str());
    else       std::printf("* SYMMETRISABLE   symmetric representative written to %s\n",
                           out.c_str());
    return 0;
}
