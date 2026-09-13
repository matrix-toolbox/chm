// =============================================================================
// chm.hpp -- Complex Hadamard Matrices: shared core library
//
// C++ port of the Octave toolbox (ud.m, PD*.m, sinkhorn.m, ss.m, sh.m,
// dephase.m, isBH.m, LF.m, SL.m, SL3.m, reshuffle.m, Tx.m, nh.m, getUnique.m,
// summary.m).
//
// 2026-09-04  W. Bruzda  <w.bruzda@cft.edu.pl>
//
// Conventions used throughout
// ---------------------------
//   * A complex Hadamard matrix (CHM) of order N satisfies
//         |H_jk| = 1           for all j, k          -> n1(H) = 0
//         H H^dagger = N I_N                         -> nh(H) = 0
//   * "core phases" A is the (N-1)x(N-1) real array such that the dephased
//     matrix reads
//         H = [ 1  1 ... 1 ]
//             [ 1             ]
//             [ :  exp(2 pi i A) ]
//     with phases normalised to the unit interval, A_jk in [0, 1).
//   * Matrices are stored ROW-MAJOR.  LAPACKE is called with
//     LAPACK_ROW_MAJOR so no transposition is ever needed.
//
// Requires: LAPACKE (Ubuntu: liblapacke-dev libopenblas-dev)
// See COMPILE.txt for the exact build lines.
// =============================================================================

#ifndef CHM_HPP
#define CHM_HPP

// -----------------------------------------------------------------------------
// Version handshake -- guards against a MIXED source tree
// -----------------------------------------------------------------------------
// The shared library gained exe_dir / run_int / name_after_matrix on 2026-09-08;
// before that those helpers lived as file-static functions inside get_chm.cpp.
// Compiling an old get_chm.cpp against this header therefore makes both visible
// (get_chm.cpp says "using namespace chm;") and every call to them is ambiguous,
// which g++ reports as a wall of overload-resolution notes that say nothing
// about the actual cause.
//
// So each program states the version it was written for, and the two sides
// check each other:
//
//   * every .cpp does  #define CHM_PROGRAM_ABI <date>  BEFORE including this
//     header, and the test below rejects a .cpp that is older than the header;
//   * every .cpp also tests CHM_ABI after the include, which rejects a header
//     older than the .cpp.
//
// Either way the message names the problem: the files do not come from the same
// archive.  Extract one archive into one empty directory.
#define CHM_ABI 20260908

#if !defined(CHM_PROGRAM_ABI)
#error "MIXED SOURCE TREE: this .cpp predates chm.hpp (it never declared CHM_PROGRAM_ABI). Delete the directory and extract one archive into it; see COMPILE.txt."
#elif CHM_PROGRAM_ABI < CHM_ABI
#error "MIXED SOURCE TREE: this .cpp is older than chm.hpp. Delete the directory and extract one archive into it; see COMPILE.txt."
#endif

#include <complex>
#include <cstdint>
#include <string>
#include <vector>

namespace chm {

using cd = std::complex<double>;

constexpr double PI  = 3.14159265358979323846;
constexpr double TAU = 6.28318530717958647692;  // 2*pi

// -----------------------------------------------------------------------------
// Square complex matrix, row-major.
// -----------------------------------------------------------------------------
struct CMat {
    int n = 0;
    std::vector<cd> a;

    CMat() = default;
    explicit CMat(int n_, cd fill = cd(0.0, 0.0)) : n(n_), a(std::size_t(n_) * n_, fill) {}

    inline cd&       operator()(int i, int j)       { return a[std::size_t(i) * n + j]; }
    inline const cd& operator()(int i, int j) const { return a[std::size_t(i) * n + j]; }

    static CMat identity(int n);
    static CMat fourier(int n);            // F_jk = exp(2 pi i j k / n), unnormalised
};

// General real matrix, row-major (used for the defect's system matrix R).
struct RMat {
    int rows = 0, cols = 0;
    std::vector<double> a;

    RMat() = default;
    RMat(int r, int c) : rows(r), cols(c), a(std::size_t(r) * c, 0.0) {}

    inline double&       operator()(int i, int j)       { return a[std::size_t(i) * cols + j]; }
    inline const double& operator()(int i, int j) const { return a[std::size_t(i) * cols + j]; }
};

// -----------------------------------------------------------------------------
// Norms and elementary transformations
// -----------------------------------------------------------------------------
double n1(const CMat& X);                  // || |X| - 1 ||_F        (unimodularity)
double nh(const CMat& X);                  // || X X^dag - N I ||_F  (orthogonality)

CMat   transpose(const CMat& X);           // X^T
CMat   adjoint(const CMat& X);             // X^dag
CMat   dephase(const CMat& X);             // first row and column set to ones
CMat   reshuffle(const CMat& X);           // X^R   -- requires N = d*d
CMat   partial_transpose(const CMat& X, int sys);  // X^Gamma, sys in {1,2}; N = d*d

double frob(const CMat& X);
double sym_error(const CMat& X);           // || X - X^T ||_F
double herm_error(const CMat& X);          // || X - X^dag ||_F

// Extract core phases in [0,1) from a dephased matrix; returns (N-1)x(N-1).
std::vector<std::vector<double>> core_phases(const CMat& H);

// Rebuild the full NxN dephased matrix from an (N-1)x(N-1) core-phase array.
CMat from_core_phases(const std::vector<std::vector<double>>& A);

// Build an NxN matrix from a full NxN phase array: H = exp(2 pi i A).
CMat from_full_phases(const std::vector<std::vector<double>>& A);

// -----------------------------------------------------------------------------
// Linear algebra (LAPACKE backend)
// -----------------------------------------------------------------------------
// Singular values only (LAPACK ?gesdd with jobz = 'N').
std::vector<double> singular_values(const RMat& A);
std::vector<double> singular_values(const CMat& A);

// Unitary polar factor U of X = U P, computed as U = W V^dag from the SVD.
// Numerically the stable replacement for PD.m / PD_SVD.m / PD_SVDi.m.
CMat polar_unitary(const CMat& X);

// Numerical rank from singular values.
//   tol > 0  : count of s_i > tol                     (Octave ud(..,'S',tol))
//   tol <= 0 : count of s_i > max(m,n)*eps*s_max      (Octave rank(), method 'R')
int rank_from_sv(const std::vector<double>& s, int rows, int cols, double tol);

// -----------------------------------------------------------------------------
// Generation of CHM
// -----------------------------------------------------------------------------
enum class Method   { Sinkhorn, RWCP };
enum class Symmetry { None, Symmetric, Hermitian };

struct GenOptions {
    int      size       = 7;
    Method   method     = Method::Sinkhorn;
    Symmetry symmetry   = Symmetry::None;
    long long iterations = 10000;   // max iterations per restart (--iterations)
    long long restarts  = 0;        // 0 = unlimited restarts until success
    std::uint64_t seed  = 0;        // 0 = draw a nondeterministic seed
    double   eps        = 1e-13;    // convergence threshold on n1 and nh
    int      threads    = 0;        // 0 = OpenMP default; used by rwcp restarts
    double   max_seconds = 0.0;     // 0 = no wall-clock limit
    int      kicks      = 8;        // perturbations tried before a full restart
    bool     polish     = true;     // sharpen a near-miss Sinkhorn iterate
    bool     verbose    = false;
};

struct GenResult {
    CMat      H;
    bool      ok        = false;
    long long iterations = 0;       // iterations spent in the successful restart
    long long restarts  = 0;        // number of restarts performed
    double    n1_val    = 0.0;
    double    nh_val    = 0.0;
    std::uint64_t seed  = 0;        // seed actually used (for reproducibility)
};

GenResult generate(const GenOptions& opt);

// A Hermitian CHM of order N has H^2 = N I, hence eigenvalues +/- sqrt(N), hence
//     tr H = sqrt(N) * (n_+ - n_-) ,
// while tr H = sum_j H_jj is a sum of N real unimodular numbers, i.e. an integer.
// For N not a perfect square sqrt(N) is irrational, which forces n_+ = n_- and
// therefore N even.  So Hermitian CHM can only exist when N is even or a
// perfect square -- there is nothing to search for at N = 3, 5, 7, 11, ...
bool hermitian_chm_possible(int n);

// -----------------------------------------------------------------------------
// Analysis
// -----------------------------------------------------------------------------
// Dephased defect d(U) = (N-1)^2 - rank(R).  Port of ud.m.
//   method 'R' : exact rank of R          (tol ignored)
//   method 'S' : rank = #{ sv(R) > tol }  (recommended for numerical matrices)
//   method 'T' : tangent-space method
int defect(const CMat& U, char method = 'S', double tol = 1e-8);

// Haagerup invariants  Lambda(H) = { H_ij H_kl conj(H_il) conj(H_kj) }.
// Returns the unique representatives, sorted by argument then modulus.
std::vector<cd> haagerup_set(const CMat& H, double eps = 1e-8);

// Smallest q, 1 < q <= qmax, with H_jk^q = 1 for all entries; 0 if none.
// Port of isBH.m, but tested through phases (no accumulation of round-off).
int butson_q(const CMat& H, int qmax = 10000, double tol = 1e-8);

// LOG-form: L_jk = round(q * arg(H_jk) / 2pi) mod q.  Port of LF.m.
std::vector<std::vector<int>> log_form(const CMat& H, int q);

// Linear ("singular") entropy, rescaled so that SL(unitary) = 1.  Port of SL.m.
double linear_entropy(const CMat& M);

struct SL3 {
    bool   defined   = false;      // false unless N is a perfect square
    int    d         = 0;          // N = d*d
    double s_H = 0.0, s_R = 0.0, s_G = 0.0;
    bool   r_dual = false, g_dual = false, two_unitary = false;
};
SL3 linear_entropy_triplet(const CMat& H, double tol = 1e-10);

// -----------------------------------------------------------------------------
// I/O
// -----------------------------------------------------------------------------
std::string timestamp();                         // local time, YYYYMMDDThhmmss

// The timestamp has one-second granularity, so two runs started inside the same
// second would collide.  Returns `path` when it is free, otherwise the first of
// path_1, path_2, ... that is (extension kept where there is one).
std::string free_path(const std::string& path);

enum class Layout { Core, Full, Auto };

struct LoadResult {
    CMat   H;
    int    array_size = 0;   // side length of the array found in the file
    Layout used       = Layout::Core;
    bool   ok         = false;
    std::string error;
};

// Read a real phase array from an Octave m-file (or a bare numeric .data file)
// and turn it into the complex matrix H.
//   Layout::Core -> array is the (N-1)x(N-1) core, H is bordered with ones
//   Layout::Full -> array is the full NxN phase matrix
//   Layout::Auto -> try Core, fall back to Full if that is closer to a CHM
LoadResult load_phase_file(const std::string& path, Layout layout = Layout::Core);

// Write the dephased core phases of H as a callable Octave function whose
// name equals the file's base name.
// digits = significant digits per phase; 17 round-trips a double exactly.
bool write_core_m(const std::string& path, const CMat& H, const std::string& stamp,
                  const std::string& provenance, int digits = 17);

// Write the unique Haagerup invariants, one "re<TAB>im" pair per line.
bool write_lambda_data(const std::string& path, const std::vector<cd>& lambda);

// -----------------------------------------------------------------------------
// Naming a result file after the matrix it holds
// -----------------------------------------------------------------------------
// Directory holding the running binary, with a trailing '/'.  The sibling
// programs ./defect and ./lambda are looked for there, so the toolbox works
// from any working directory -- and also when the binary was found on PATH, in
// which case argv[0] carries no directory at all and /proc/self/exe answers.
std::string exe_dir(const std::string& argv0);

// Shell-quote a path so that spaces and quotes survive popen().
std::string shell_quote(const std::string& s);

// Run `cmd` and expect a single integer on its stdout.  False if the program is
// missing, exits non-zero, or prints anything else.
bool run_int(const std::string& cmd, long long& out);

// S if H is symmetric, H if Hermitian, Y if neither (tolerance on the Frobenius
// distance).  Symmetric wins when a matrix is somehow both.
char symmetry_letter(const CMat& H, double tol = 1e-7);

// Rename `path` (already written, holding H) to  <x>H_<N>_<d>_<L>_<stamp>.m by
// re-writing it under the new name -- Octave resolves a function by the file's
// base name, so a plain rename would leave the file uncallable.
//
// d and L come from running ./defect and ./lambda in `bindir`, unless they are
// passed in already (>= 0), which is what a program that knows them from an
// input file name does: both are invariants of the monomial equivalence class.
// Progress and failures are reported on `log` (nullptr = silent).  If either
// number cannot be had, `path` is left untouched and returned unchanged.
std::string name_after_matrix(const std::string& path, const CMat& H,
                              const std::string& stamp, const std::string& provenance,
                              int digits, const std::string& bindir,
                              std::FILE* log, long long d = -1, long long L = -1);

// Parse  <x>H_<N>_<d>_<L>_<stamp>.m , the name written by get_chm.  Returns
// false unless every field is present and well formed; on success N, d and L
// are filled in.  Used to reuse the two invariants instead of recomputing them.
bool parse_matrix_name(const std::string& path, int& N, long long& d, long long& L);

}  // namespace chm

#endif  // CHM_HPP
