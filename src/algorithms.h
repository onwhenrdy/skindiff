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

    // Fused Crank-Nicolson sub-step:  vec <- lhs^{-1} * rhs * vec
    // Equivalent to inlineMultiply(rhs, vec) followed by thomasReUseIP(lhs, vec),
    // but does one pass instead of two: the multiplied value stays in a register
    // and never round-trips through the cache between the multiply and the
    // forward solve. The backward sweep is unchanged (it can't fuse because it
    // depends on the fully-forwarded vector).
    //
    // lhs is mutated into the prepared form on first call (same convention as
    // thomasReUseIP). rhs is read-only.
    inline void crankNicolsonStepIP(const TDMatrix& rhs_mat, TDMatrix& lhs,
                                    std::vector<double>& vec)
    {
        const auto size = lhs.size();
        assert(size > 1);
        assert(static_cast<std::size_t>(size) == vec.size());
        assert(rhs_mat.size() == size);

        auto& c_star  = lhs.fullUpper();
        auto& c_diag  = lhs.fullDiag();
        auto& c_lower = lhs.fullLower();

        if (!lhs.isPrepared())
        {
            c_star[0] = c_star[0] / c_diag[0];
            for (int i = 1; i < size - 1; ++i)
            {
                c_diag[i] = c_diag[i] - c_star[i - 1] * c_lower[i - 1];
                c_star[i] = c_star[i] / c_diag[i];
            }
            c_diag[size - 1] = c_diag[size - 1] - c_star[size - 2] * c_lower[size - 2];
            for (int i = 0; i < size; ++i) c_diag[i] = 1.0 / c_diag[i];
            lhs.setPrepared(true);
        }

        const auto& m_diag  = rhs_mat.fullDiag();
        const auto& m_upper = rhs_mat.fullUpper();
        const auto& m_lower = rhs_mat.fullLower();

        // Boundary cell 0: M*vec uses only diag and upper; forward solve at row
        // 0 has no sub-diagonal contribution.
        double tmp_prev      = vec[0];                    // old vec[i-1] for the M*vec recurrence
        const double mul_0   = m_diag[0] * vec[0] + m_upper[0] * vec[1];
        vec[0]               = mul_0 * c_diag[0];

        // Interior: compute M*vec at index i and apply the forward solve in
        // the same step. tmp_prev carries the unmodified vec[i-1] forward
        // through the M*vec recurrence (same trick as inlineMultiply).
        for (int i = 1; i < size - 1; ++i)
        {
            const double old_vec_i = vec[i];
            const double mul_i = m_lower[i - 1] * tmp_prev
                               + m_diag[i]      * old_vec_i
                               + m_upper[i]     * vec[i + 1];
            vec[i] = (mul_i - vec[i - 1] * c_lower[i - 1]) * c_diag[i];
            tmp_prev = old_vec_i;
        }

        // Last cell: M*vec uses only lower and diag.
        {
            const auto idx = size - 1;
            const double mul_last = m_lower[idx - 1] * tmp_prev + m_diag[idx] * vec[idx];
            vec[idx] = (mul_last - vec[idx - 1] * c_lower[idx - 1]) * c_diag[idx];
        }

        // Backward sweep: unchanged.
        for (int i = size - 2; i >= 0; --i)
        {
            vec[i] = vec[i] - c_star[i] * vec[i + 1];
        }
    }
}

#endif  // SC_ALGORITHMS_H
