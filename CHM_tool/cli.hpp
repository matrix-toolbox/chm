// =============================================================================
// cli.hpp -- minimal command-line helpers shared by the four programs
//
// 2026-09-04  W. Bruzda  <w.bruzda@cft.edu.pl>
//
// Header-only, so it never appears on a compile line of its own.
// =============================================================================

#ifndef CHM_CLI_HPP
#define CHM_CLI_HPP

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace cli {

struct Args {
    std::string prog;
    std::vector<std::string> tok;
    std::string error;

    bool has(const std::string& flag) const {
        for (const auto& t : tok) if (t == flag) return true;
        return false;
    }

    // Value of "--flag value"; returns false when the flag is absent.
    bool value(const std::string& flag, std::string& out) {
        for (std::size_t i = 0; i < tok.size(); ++i) {
            if (tok[i] != flag) continue;
            if (i + 1 >= tok.size() || tok[i + 1].rfind("--", 0) == 0) {
                error = "flag '" + flag + "' needs a value";
                return false;
            }
            out = tok[i + 1];
            return true;
        }
        return false;
    }

    // Accepts plain integers as well as the "1e+4" form used in the spec.
    bool integer(const std::string& flag, long long& out) {
        std::string s;
        if (!value(flag, s)) return false;
        char* end = nullptr;
        const double v = std::strtod(s.c_str(), &end);
        if (end == s.c_str() || (end && *end != '\0') || !std::isfinite(v)) {
            error = "flag '" + flag + "': '" + s + "' is not a number";
            return false;
        }
        out = static_cast<long long>(std::llround(v));
        return true;
    }

    bool real(const std::string& flag, double& out) {
        std::string s;
        if (!value(flag, s)) return false;
        char* end = nullptr;
        const double v = std::strtod(s.c_str(), &end);
        if (end == s.c_str() || (end && *end != '\0')) {
            error = "flag '" + flag + "': '" + s + "' is not a number";
            return false;
        }
        out = v;
        return true;
    }

    // Reject anything that looks like a flag but is not in the known list, so
    // that a typo such as "--symetry" fails loudly instead of being ignored.
    bool check_known(const std::vector<std::string>& known,
                     const std::vector<std::string>& valued) {
        for (std::size_t i = 0; i < tok.size(); ++i) {
            if (tok[i].rfind("--", 0) != 0 && tok[i].rfind("-", 0) != 0) continue;
            bool ok = false;
            for (const auto& k : known) if (tok[i] == k) { ok = true; break; }
            if (!ok) { error = "unknown flag '" + tok[i] + "'"; return false; }
            for (const auto& k : valued) if (tok[i] == k) { ++i; break; }
        }
        return true;
    }
};

inline Args parse(int argc, char** argv) {
    Args a;
    a.prog = (argc > 0 && argv[0]) ? argv[0] : "program";
    for (int i = 1; i < argc; ++i) a.tok.emplace_back(argv[i]);
    return a;
}

}  // namespace cli

#endif  // CHM_CLI_HPP
