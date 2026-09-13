// =============================================================================
// summary.cpp -- check if a given matrix is CHM
//
// 2026-09-04  W. Bruzda  <w.bruzda@cft.edu.pl>
//
// Reads the phase array produced by get_chm, rebuilds H = exp(2 pi i A) and
// reports:  the CHM test, the Butson class, symmetry, the cardinality of the
// Haagerup invariants, the defect and the triplet of linear entropies.
// C++ port of summary.m.
//
// Build:
//   g++ -std=c++17 -O3 -march=znver3 -mtune=znver3 -fopenmp -flto
//       -o summary summary.cpp chm.cpp -llapacke -llapack -lopenblas -lm
// =============================================================================

#define CHM_PROGRAM_ABI 20260908   // must match CHM_ABI in chm.hpp
#include "chm.hpp"
#if CHM_ABI < CHM_PROGRAM_ABI
#error "MIXED SOURCE TREE: chm.hpp is older than this .cpp. Delete the directory and extract one archive into it; see COMPILE.txt."
#endif

#include "cli.hpp"

#include <cstdio>
#include <string>

using namespace chm;

static void usage(const char* prog) {
    std::printf(
        "usage: %s --file FILE [--symmetry] [--defect] [--butson] [--lambda]\n"
        "          [--SL3] [--all] [--layout core|full|auto] [--eps E]\n"
        "          [--tol T] [--lambda-eps E] [--qmax Q] [--butson-tol T] [--log-form]\n"
        "\n"
        "  --file FILE     phase array written by get_chm (mandatory)\n"
        "  --symmetry      test H = H^T and H = H^dag        (always done; default)\n"
        "  --defect        value of the dephased defect\n"
        "  --butson        smallest q with H_jk^q = 1 for all entries\n"
        "  --lambda        cardinality of the Haagerup invariants\n"
        "  --SL3           linear entropies of H, H^R and H^G; needs N = d*d\n"
        "  --all           everything above\n"
        "  --layout L      core | full | auto                       (default core)\n"
        "  --eps E         CHM threshold on nh and n1               (default 1e-7)\n"
        "  --tol T         singular-value threshold for the defect  (default 1e-8)\n"
        "  --lambda-eps E  uniqueness threshold for the invariants   (default 1e-8)\n"
        "  --qmax Q        largest Butson exponent tested           (default 10000)\n"
        "  --butson-tol T  tolerance of the Butson test             (default 1e-8)\n"
        "  --log-form      print the LOG-form matrix when H is of Butson type\n"
        "\n"
        "flags may be given in any order\n"
        "\n"
        "exit: 0 = the matrix is a CHM, 1 = it is not, 2 = bad usage/unreadable file\n",
        prog);
}

int main(int argc, char** argv) {
    cli::Args args = cli::parse(argc, argv);

    const std::vector<std::string> known = {
        "--file",   "--symmetry",    "--defect",   "--butson", "--lambda", "--SL3",
        "--all",    "--layout",      "--eps",      "--tol",    "--lambda-eps",
        "--qmax",   "--butson-tol",  "--log-form", "--help",   "-h"};
    const std::vector<std::string> valued = {"--file", "--layout",     "--eps",
                                             "--tol",  "--lambda-eps", "--qmax",
                                             "--butson-tol"};

    if (args.has("--help") || args.has("-h")) { usage(args.prog.c_str()); return 0; }

    // "$ ./summary" -- no arguments prints a warning
    if (args.tok.empty()) {
        std::fprintf(stderr, "* warning: no arguments given; --file is mandatory\n");
        usage(args.prog.c_str());
        return 2;
    }
    if (!args.check_known(known, valued)) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), args.error.c_str());
        return 2;
    }

    std::string file, s;
    if (!args.value("--file", file)) {
        std::fprintf(stderr, "* warning: --file is mandatory\n");
        return 2;
    }

    const bool all       = args.has("--all");
    const bool do_defect = all || args.has("--defect");
    const bool do_butson = all || args.has("--butson");
    const bool do_lambda = all || args.has("--lambda");
    const bool do_sl3    = all || args.has("--SL3");
    const bool do_log    = args.has("--log-form");

    Layout layout = Layout::Core;
    if (args.value("--layout", s)) {
        if      (s == "core") layout = Layout::Core;
        else if (s == "full") layout = Layout::Full;
        else if (s == "auto") layout = Layout::Auto;
        else { std::fprintf(stderr, "%s: unknown layout '%s'\n", args.prog.c_str(), s.c_str()); return 2; }
    }

    double    eps = 1e-7, tol = 1e-8, leps = 1e-8, btol = 1e-8;
    long long qmax = 10000;
    args.real("--eps", eps);
    args.real("--tol", tol);
    args.real("--lambda-eps", leps);
    args.real("--butson-tol", btol);
    args.integer("--qmax", qmax);
    if (!args.error.empty()) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), args.error.c_str());
        return 2;
    }

    const LoadResult L = load_phase_file(file, layout);
    if (!L.ok) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), L.error.c_str());
        return 2;
    }
    const CMat& H = L.H;
    const int   N = H.n;

    const double nf = nh(H);
    const double v1 = n1(H);

    // ---- the CHM test gates everything else --------------------------------
    if (nf > eps || v1 > eps) {
        std::printf("* this is not a CHM!\n");
        std::printf("* N = %d, nf = %g, n1 = %g (threshold %g)\n", N, nf, v1, eps);
        return 1;
    }
    std::printf("* this is a CHM (nf = %g, n1 = %g)\n", nf, v1);

    // ---- Butson class -------------------------------------------------------
    if (do_butson) {
        const int q = butson_q(H, int(qmax), btol);
        if (q) {
            std::printf("* this is a BH(%d, %d)\n", N, q);
            if (do_log) {
                const auto Lg = log_form(H, q);
                for (int i = 0; i < N; ++i) {
                    std::printf("  ");
                    for (int j = 0; j < N; ++j) std::printf(" %*d", (q > 9 ? 3 : 2), Lg[std::size_t(i)][std::size_t(j)]);
                    std::printf("\n");
                }
            }
        } else {
            std::printf("* this is probably not a BH-matrix for any q < %lld\n", qmax + 1);
        }
    }

    // ---- symmetry (default, always reported) --------------------------------
    const double se = sym_error(H);
    const double he = herm_error(H);
    if (se <= eps) std::printf("* this is a symmetric matrix\n");
    if (he <= eps) std::printf("* this is a Hermitian matrix\n");

    // ---- Haagerup invariants ------------------------------------------------
    if (do_lambda) {
        const std::vector<cd> lam = haagerup_set(H, leps);
        std::printf("* #L = %zu\n", lam.size());
    }

    // ---- defect -------------------------------------------------------------
    if (do_defect) {
        std::printf("* ud = %d\n", defect(H, 'S', tol));
    }

    // ---- linear entropies ---------------------------------------------------
    if (do_sl3) {
        const SL3 e = linear_entropy_triplet(H);
        if (!e.defined) {
            std::printf("* SL(H) = n/a  (H^R and H^G need N = d*d; here N = %d)\n", N);
        } else {
            std::printf("* SL(H) = %.12g %.12g %.12g\n", e.s_H, e.s_R, e.s_G);
            if (e.two_unitary)   std::printf("* 2-unitary matrix!\n");
            else if (e.r_dual)   std::printf("* R-dual matrix!\n");
            else if (e.g_dual)   std::printf("* G-dual matrix!\n");
        }
    }

    return 0;
}
