// =============================================================================
// lambda.cpp -- cardinality of the (unique) Haagerup invariants
//
// 2026-09-04  W. Bruzda  <w.bruzda@cft.edu.pl>
//
//   Lambda(H) = { H_ij H_kl conj(H_il) conj(H_kj) : i,j,k,l = 1..N }
//
// The unique set is written to lambda_YYYYMMDDThhmmss.data, one invariant per
// line as two real numbers (real part, imaginary part) separated by a tab.
//
// Build:
//   g++ -std=c++17 -O3 -march=znver3 -mtune=znver3 -fopenmp -flto
//       -o lambda lambda.cpp chm.cpp -llapacke -llapack -lopenblas -lm
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
        "usage: %s --file FILE [--layout core|full|auto] [--eps E]\n"
        "          [--out FILE] [--no-write] [--bare]\n"
        "\n"
        "  --file FILE     phase array written by get_chm (mandatory)\n"
        "  --layout L      core | full | auto                     (default core)\n"
        "  --eps E         two invariants count as equal when their distance\n"
        "                  is below E                             (default 1e-8)\n"
        "  --out FILE      output path    (default lambda_YYYYMMDDThhmmss.data)\n"
        "  --no-write      report the cardinality without writing the set\n"
        "  --bare          print the cardinality only, nothing else\n"
        "\n"
        "exit: 0 = value printed, 2 = bad usage or unreadable file\n",
        prog);
}

int main(int argc, char** argv) {
    cli::Args args = cli::parse(argc, argv);

    const std::vector<std::string> known  = {"--file", "--layout",   "--eps",  "--out",
                                             "--no-write", "--bare", "--help", "-h"};
    const std::vector<std::string> valued = {"--file", "--layout", "--eps", "--out"};

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

    std::string file, out, s;
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

    double eps = 1e-8;
    args.real("--eps", eps);
    args.value("--out", out);
    if (!args.error.empty()) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), args.error.c_str());
        return 2;
    }
    if (!(eps > 0.0)) {
        std::fprintf(stderr, "%s: --eps must be positive\n", args.prog.c_str());
        return 2;
    }

    const LoadResult L = load_phase_file(file, layout);
    if (!L.ok) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), L.error.c_str());
        return 2;
    }

    const std::vector<cd> lam = haagerup_set(L.H, eps);

    const bool bare  = args.has("--bare");
    const bool write = !args.has("--no-write");

    if (write) {
        if (out.empty()) out = free_path("lambda_" + timestamp() + ".data");
        if (!write_lambda_data(out, lam)) {
            std::fprintf(stderr, "%s: cannot write '%s'\n", args.prog.c_str(), out.c_str());
            return 2;
        }
    }

    if (bare) {
        std::printf("%zu\n", lam.size());
    } else {
        std::printf("* N  = %d (%s layout)\n", L.H.n, L.used == Layout::Core ? "core" : "full");
        std::printf("* nh = %.6e, n1 = %.6e\n", nh(L.H), n1(L.H));
        std::printf("* eps = %g\n", eps);
        std::printf("* #L = %zu\n", lam.size());
        if (write) std::printf("* unique invariants written to %s\n", out.c_str());
    }
    return 0;
}
