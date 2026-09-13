// =============================================================================
// defect.cpp -- calculate the defect of a unitary matrix
//
// 2026-09-04  W. Bruzda  <w.bruzda@cft.edu.pl>
//
// The dephased defect  d(U) = (N-1)^2 - rank(R),  where R is the real system
// matrix of the linearised phase-perturbation conditions.  See
//   W. Tadej, K. Zyczkowski, Linear Algebra Appl. 429, 447-481 (2008).
// C++ port of ud.m.
//
// Build:
//   g++ -std=c++17 -O3 -march=znver3 -mtune=znver3 -fopenmp -flto
//       -o defect defect.cpp chm.cpp -llapacke -llapack -lopenblas -lm
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
        "usage: %s --file FILE [--layout core|full|auto] [--method R|S|T]\n"
        "          [--tol T] [--bare]\n"
        "\n"
        "  --file FILE     phase array written by get_chm (mandatory)\n"
        "  --layout L      core = FILE holds the (N-1)x(N-1) core   (default core)\n"
        "                  full = FILE holds the full NxN phases\n"
        "                  auto = pick whichever is closer to a CHM\n"
        "  --method M      R = exact rank of R\n"
        "                  S = rank as #{ singular values > --tol }  (default)\n"
        "                  T = tangent-space method; R and T are the same matrix\n"
        "                      up to transposition, so they always agree\n"
        "  --tol T         singular-value threshold for method S     (default 1e-8)\n"
        "  --bare          print the integer only, nothing else\n"
        "\n"
        "exit: 0 = value printed, 2 = bad usage or unreadable file\n",
        prog);
}

int main(int argc, char** argv) {
    cli::Args args = cli::parse(argc, argv);

    const std::vector<std::string> known  = {"--file", "--layout", "--method", "--tol",
                                             "--bare", "--help",   "-h"};
    const std::vector<std::string> valued = {"--file", "--layout", "--method", "--tol"};

    if (args.has("--help") || args.has("-h")) { usage(args.prog.c_str()); return 0; }
    if (args.tok.empty()) {
        std::fprintf(stderr, "* warning: --file is mandatory\n");
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

    Layout layout = Layout::Core;
    if (args.value("--layout", s)) {
        if      (s == "core") layout = Layout::Core;
        else if (s == "full") layout = Layout::Full;
        else if (s == "auto") layout = Layout::Auto;
        else { std::fprintf(stderr, "%s: unknown layout '%s'\n", args.prog.c_str(), s.c_str()); return 2; }
    }

    char method = 'S';
    if (args.value("--method", s)) {
        if (s.size() != 1 || (s[0] != 'R' && s[0] != 'S' && s[0] != 'T')) {
            std::fprintf(stderr, "%s: unknown method '%s' (use R, S or T)\n",
                         args.prog.c_str(), s.c_str());
            return 2;
        }
        method = s[0];
    }

    double tol = 1e-8;
    args.real("--tol", tol);
    if (!args.error.empty()) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), args.error.c_str());
        return 2;
    }

    const LoadResult L = load_phase_file(file, layout);
    if (!L.ok) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), L.error.c_str());
        return 2;
    }

    const int d = defect(L.H, method, tol);

    if (args.has("--bare")) {
        std::printf("%d\n", d);
    } else {
        std::printf("* N  = %d (%s layout, %dx%d array in file)\n", L.H.n,
                    L.used == Layout::Core ? "core" : "full", L.array_size, L.array_size);
        std::printf("* nh = %.6e, n1 = %.6e\n", nh(L.H), n1(L.H));
        if (method == 'S') std::printf("* method = S, tol = %g\n", tol);
        else               std::printf("* method = %c\n", method);
        std::printf("* ud = %d\n", d);
    }
    return 0;
}
