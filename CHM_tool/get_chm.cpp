// =============================================================================
// get_chm.cpp -- return a complex Hadamard matrix of given size
//
// 2026-09-04  W. Bruzda  <w.bruzda@cft.edu.pl>
//
// Writes the dephased matrix as its (N-1)x(N-1) array of core phases into
//
//     <x>H_<N>_<d>_<L>_<YYYYMMDDThhmmss>.m
//
// a callable Octave function, where
//     x = S  the matrix is symmetric      (H = H^T)
//     x = H  the matrix is Hermitian      (H = H^dag, and not symmetric)
//     x = Y  neither
//     N      the order of the matrix
//     d      the dephased defect, obtained by running  ./defect
//     L      the Haagerup cardinality,  obtained by running  ./lambda
//
// If either sibling program cannot be run the file keeps the neutral name
// A_<YYYYMMDDThhmmss>.m .  An explicit --out is never renamed.
//
// Build:
//   g++ -std=c++17 -O3 -march=znver3 -mtune=znver3 -fopenmp -flto
//       -o get_chm get_chm.cpp chm.cpp -llapacke -llapack -lopenblas -lm
// =============================================================================

#define CHM_PROGRAM_ABI 20260908   // must match CHM_ABI in chm.hpp
#include "chm.hpp"
#if CHM_ABI < CHM_PROGRAM_ABI
#error "MIXED SOURCE TREE: chm.hpp is older than this .cpp. Delete the directory and extract one archive into it; see COMPILE.txt."
#endif

#include "cli.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

using namespace chm;

// -----------------------------------------------------------------------------
// naming the output file:   <x>H_<N>_<d>_<L>_<YYYYMMDDThhmmss>.m
//
//   x = S  symmetric,  H  Hermitian,  Y  neither
//   d = the dephased defect, from ./defect
//   L = the cardinality of the Haagerup invariants, from ./lambda
//
// The symmetry letter is free -- it is already known from the matrix -- but d
// and L are obtained by running the two sibling programs, as specified.  If
// either cannot be run the file keeps the neutral name A_<stamp>.m.
//
// The machinery lives in chm.cpp (exe_dir / run_int / name_after_matrix), since
// symmetrize writes result files the same way.
// -----------------------------------------------------------------------------

static void usage(const char* prog) {
    std::printf(
        "usage: %s [--size N] [--method sinkhorn|rwcp] [--iterations M]\n"
        "          [--symmetry s|h|none] [--seed S] [--restarts R] [--eps E]\n"
        "          [--digits D] [--threads T] [--out FILE] [--quiet] [--verbose]\n"
        "\n"
        "  --size N        order of the matrix                    (default 7)\n"
        "  --method M      sinkhorn | rwcp                        (default sinkhorn)\n"
        "                    sinkhorn = alternating projection between the\n"
        "                               unimodular and the unitary manifold\n"
        "                    rwcp     = random walk over core phases\n"
        "  --iterations M  max iterations per restart, 1e+4 accepted (default 1e+4)\n"
        "                  --iteration is accepted as an alias\n"
        "  --symmetry X    s = symmetric, h = Hermitian, none      (default none)\n"
        "  --seed S        PRNG seed; 0 draws one from the OS      (default 0)\n"
        "  --restarts R    max restarts, 0 = until success         (default 0)\n"
        "  --eps E         convergence threshold on n1 and nh      (default 1e-13)\n"
        "  --digits D      significant digits written to the file  (default 17)\n"
        "  --threads T     restart threads, 0 = all cores          (default 0)\n"
        "  --kicks K       phase perturbations tried before giving up on a\n"
        "                  restart, 0 = restart immediately        (default 8)\n"
        "  --max-seconds S wall-clock budget, 0 = unlimited        (default 0)\n"
        "  --no-polish     skip the exact-coordinate-descent stage that\n"
        "                  sharpens a near-miss sinkhorn iterate\n"
        "  --out FILE      output path; suppresses the automatic naming\n"
        "                  (default <x>H_<N>_<d>_<L>_<YYYYMMDDThhmmss>.m, with\n"
        "                   x = S/H/Y for symmetric/Hermitian/neither, and d, L\n"
        "                   taken from ./defect and ./lambda; if either program\n"
        "                   is missing the name stays A_<stamp>.m)\n"
        "  --quiet         print the output file name only\n"
        "  --verbose       one dot per failed restart, on stderr\n"
        "  --help, -h      this message\n"
        "\n"
        "exit: 0 = a CHM was found and written, 1 = search failed, 2 = bad usage\n",
        prog);
}

int main(int argc, char** argv) {
    cli::Args args = cli::parse(argc, argv);

    const std::vector<std::string> known = {
        "--size", "--method", "--iterations", "--iteration", "--symmetry", "--seed",
        "--restarts", "--eps", "--digits", "--threads", "--out", "--quiet",
        "--verbose", "--kicks", "--max-seconds", "--no-polish", "--help", "-h"};
    const std::vector<std::string> valued = {
        "--size", "--method", "--iterations", "--iteration", "--symmetry", "--seed",
        "--restarts", "--eps", "--digits", "--threads", "--out", "--kicks",
        "--max-seconds"};

    if (args.has("--help") || args.has("-h")) { usage(args.prog.c_str()); return 0; }
    if (!args.check_known(known, valued)) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), args.error.c_str());
        usage(args.prog.c_str());
        return 2;
    }

    GenOptions opt;
    opt.size       = 7;         // spec default
    opt.method     = Method::Sinkhorn;
    opt.symmetry   = Symmetry::None;
    opt.iterations = 10000;     // spec default M = 1e+4
    opt.eps        = 1e-13;

    const bool quiet = args.has("--quiet");
    opt.polish       = !args.has("--no-polish");
    opt.verbose      = args.has("--verbose");
    int digits       = 17;
    std::string out;

    long long v = 0;
    double d = 0.0;
    std::string s;

    if (args.integer("--size", v))       opt.size = int(v);
    if (args.integer("--iterations", v)) opt.iterations = v;
    if (args.integer("--iteration", v))  opt.iterations = v;   // spec spells it singular
    if (args.integer("--seed", v))       opt.seed = std::uint64_t(v < 0 ? -v : v);
    if (args.integer("--restarts", v))   opt.restarts = v;
    if (args.integer("--threads", v))    opt.threads = int(v);
    if (args.integer("--digits", v))     digits = int(v);
    if (args.integer("--kicks", v))      opt.kicks = int(v);
    if (args.real("--eps", d))           opt.eps = d;
    if (args.real("--max-seconds", d))   opt.max_seconds = d;
    if (!args.error.empty()) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), args.error.c_str());
        return 2;
    }

    if (args.value("--method", s)) {
        if      (s == "sinkhorn") opt.method = Method::Sinkhorn;
        else if (s == "rwcp")     opt.method = Method::RWCP;
        else {
            std::fprintf(stderr, "%s: unknown method '%s' (use sinkhorn or rwcp)\n",
                         args.prog.c_str(), s.c_str());
            return 2;
        }
    }
    if (args.value("--symmetry", s)) {
        if      (s == "s" || s == "sym"  || s == "symmetric") opt.symmetry = Symmetry::Symmetric;
        else if (s == "h" || s == "herm" || s == "hermitian") opt.symmetry = Symmetry::Hermitian;
        else if (s == "none" || s == "n" || s == "g")         opt.symmetry = Symmetry::None;
        else {
            std::fprintf(stderr, "%s: unknown symmetry '%s' (use s, h or none)\n",
                         args.prog.c_str(), s.c_str());
            return 2;
        }
    }
    args.value("--out", out);
    if (!args.error.empty()) {
        std::fprintf(stderr, "%s: %s\n", args.prog.c_str(), args.error.c_str());
        return 2;
    }

    if (opt.size < 1) {
        std::fprintf(stderr, "%s: --size must be at least 1\n", args.prog.c_str());
        return 2;
    }
    if (opt.iterations < 1) {
        std::fprintf(stderr, "%s: --iterations must be at least 1\n", args.prog.c_str());
        return 2;
    }
    if (digits < 6 || digits > 17) {
        std::fprintf(stderr, "%s: --digits must lie between 6 and 17\n", args.prog.c_str());
        return 2;
    }
    // Refuse a search that provably cannot succeed rather than spinning for
    // ever: a Hermitian CHM has eigenvalues +/- sqrt(N) and an integer trace,
    // which is impossible unless N is even or a perfect square.
    if (opt.symmetry == Symmetry::Hermitian && !hermitian_chm_possible(opt.size)) {
        std::fprintf(stderr,
                     "%s: no Hermitian CHM of order %d exists.\n"
                     "  H = H^dag and H*H = %d*I force eigenvalues +/- sqrt(%d), so\n"
                     "  tr H = sqrt(%d)*(n+ - n-), while tr H = sum_j H_jj is an integer.\n"
                     "  For %d not a perfect square that needs n+ = n-, hence %d even.\n"
                     "  Use an even order or a perfect square.\n",
                     args.prog.c_str(), opt.size, opt.size, opt.size, opt.size, opt.size,
                     opt.size);
        return 2;
    }

    const std::string stamp = timestamp();
    // The default name starts as the neutral A_<stamp>.m and is renamed below
    // once d and L are known; an explicit --out is the caller's choice and is
    // left exactly as given.
    const bool auto_name = out.empty();
    if (auto_name) out = free_path("A_" + stamp + ".m");

    const GenResult r = generate(opt);

    if (!r.ok) {
        std::fprintf(stderr,
                     "* search failed: no CHM(%d) found within the given budget\n"
                     "  restarts = %lld, best n1 = %.3e, best nh = %.3e, seed = %llu\n",
                     opt.size, r.restarts, r.n1_val, r.nh_val,
                     static_cast<unsigned long long>(r.seed));
        if (r.n1_val > 0.0 && r.n1_val < 1e-3)
            std::fprintf(stderr,
                         "  the iterate came close: raise --iterations (the default 1e+4\n"
                         "  suffices near N = 7 but N = 16 needs ~3e+5 and N = 20 ~1e+6)\n");
        else if (opt.method == Method::Sinkhorn)
            std::fprintf(stderr,
                         "  try a larger --iterations, more --kicks, or --method rwcp\n");
        return 1;
    }

    const char* mname = (opt.method == Method::Sinkhorn) ? "sinkhorn" : "rwcp";
    const char* sname = (opt.symmetry == Symmetry::Symmetric) ? "symmetric"
                        : (opt.symmetry == Symmetry::Hermitian) ? "Hermitian" : "none";

    char prov[512];
    std::snprintf(prov, sizeof(prov),
                  "N = %d, method = %s, symmetry = %s\n"
                  "seed = %llu, iterations = %lld, restarts = %lld\n"
                  "nh = || H*H' - N*I ||_F = %.3e\n"
                  "n1 = ||  |H| - 1     ||_F = %.3e",
                  opt.size, mname, sname, static_cast<unsigned long long>(r.seed),
                  r.iterations, r.restarts, r.nh_val, r.n1_val);

    if (!write_core_m(out, r.H, stamp, prov, digits)) {
        std::fprintf(stderr, "%s: cannot write '%s'\n", args.prog.c_str(), out.c_str());
        return 1;
    }

    // ---- rename to <x>H_<N>_<d>_<L>_<stamp>.m -----------------------------------
    // Progress goes to stderr under --quiet so that stdout stays exactly one
    // line (the file name), which is what search_chm.sh and friends capture.
    std::FILE* log = quiet ? stderr : stdout;

    if (auto_name)
        out = name_after_matrix(out, r.H, stamp, prov, digits,
                                exe_dir(args.prog), log);

    if (quiet) {
        std::printf("%s\n", out.c_str());
    } else {
        std::printf("* CHM(%d) found: method = %s, symmetry = %s\n", opt.size, mname, sname);
        std::printf("* nh = %.6e, n1 = %.6e\n", r.nh_val, r.n1_val);
        std::printf("* seed = %llu, iterations = %lld, restarts = %lld\n",
                    static_cast<unsigned long long>(r.seed), r.iterations, r.restarts);
        std::printf("* || H - H^T ||_F = %.6e,  || H - H^dag ||_F = %.6e\n",
                    sym_error(r.H), herm_error(r.H));
        std::printf("* %d x %d array of core phases written to %s (%d digits)\n",
                    opt.size - 1, opt.size - 1, out.c_str(), digits);
    }
    return 0;
}
