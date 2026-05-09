#ifndef SC_ALGORITHMS_H
#define SC_ALGORITHMS_H

#include "tdmatrix.h"

#include <cassert>
#include <cstddef>
#include <vector>

namespace sc::algorithm
{
    // Solve M*x = rhs for a tri-diagonal M using the Thomas algorithm.
    // The matrix is not modified. The solution overwrites rhs.
    inline void thomasIP(const TDMatrix& matrix, std::vector<double>& rhs)
    {
        const auto size = matrix.size();
        assert(size > 0);
        assert(static_cast<std::size_t>(size) == rhs.size());

        auto c_star = matrix.fullUpper();

        c_star[0] = c_star[0] / matrix.diag(0);
        for (int i = 1; i < size - 1; ++i)
        {
            c_star[i] = c_star[i] / (matrix.diag(i) - c_star[i - 1] * matrix.lower(i - 1));
        }

        rhs[0] = rhs[0] / matrix.diag(0);
        for (int i = 1; i < size; ++i)
        {
            rhs[i] = (rhs[i] - rhs[i - 1] * matrix.lower(i - 1)) /
                     (matrix.diag(i) - c_star[i - 1] * matrix.lower(i - 1));
        }

        for (int i = size - 2; i >= 0; --i)
        {
            rhs[i] = rhs[i] - c_star[i] * rhs[i + 1];
        }
    }

    // Solve M*x = rhs reusing the LU factorization stored in M.
    // First call mutates M into the prepared form; subsequent calls with the same
    // M skip the factorization and reuse it.
    //
    // Storage convention after preparation:
    //   m_diag[i]  holds 1 / d_prepared[i]  (reciprocal -- so the forward sweep
    //                                        is multiply-only, no divides)
    //   m_upper[i] holds c_star[i]
    //   m_lower[i] holds the original sub-diagonal (unchanged)
    inline void thomasReUseIP(TDMatrix& matrix, std::vector<double>& rhs)
    {
        const auto size = matrix.size();
        assert(size > 0);
        assert(static_cast<std::size_t>(size) == rhs.size());

        auto& c_star  = matrix.fullUpper();
        auto& c_diag  = matrix.fullDiag();
        auto& c_lower = matrix.fullLower();

        if (!matrix.isPrepared())
        {
            // Fuse the c_star and c_diag updates into one pass, then invert
            // c_diag so the per-step solve uses multiplies only.
            c_star[0] = c_star[0] / c_diag[0];
            for (int i = 1; i < size - 1; ++i)
            {
                c_diag[i] = c_diag[i] - c_star[i - 1] * c_lower[i - 1];
                c_star[i] = c_star[i] / c_diag[i];
            }
            c_diag[size - 1] = c_diag[size - 1] - c_star[size - 2] * c_lower[size - 2];
            for (int i = 0; i < size; ++i) c_diag[i] = 1.0 / c_diag[i];
            matrix.setPrepared(true);
        }

        rhs[0] = rhs[0] * c_diag[0];
        for (int i = 1; i < size; ++i)
        {
            rhs[i] = (rhs[i] - rhs[i - 1] * c_lower[i - 1]) * c_diag[i];
        }

        for (int i = size - 2; i >= 0; --i)
        {
            rhs[i] = rhs[i] - c_star[i] * rhs[i + 1];
        }
    }
}

#endif  // SC_ALGORITHMS_H
