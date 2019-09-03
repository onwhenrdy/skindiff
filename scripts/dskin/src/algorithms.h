#ifndef SC_ALGORITHMS_H
#define SC_ALGORITHMS_H

#include "tdmatrix.h"
#include <cassert>
#include <cmath>

namespace sc
{
    namespace algorithm
    {
        // does modify rhs but not the matrix
        inline void thomasIP(const TDMatrix& matrix, std::vector<double>& rhs)
        {
            const auto size = matrix.size();
            assert(size > 0);
            assert(size == rhs.size());

            // copy upper band
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

            // backtracking
            for (int i = size - 2; i >= 0; --i)
            {
                rhs[i] = rhs[i] - c_star[i] * rhs[i + 1];
            }
        }

        // does modify rhs and the matrix
        inline void thomasReUseIP(TDMatrix& matrix, std::vector<double>& rhs)
        {
            const auto size = matrix.size();
            assert(size > 0);
            assert(size == rhs.size());

            auto& c_star = matrix.fullUpper();
            auto& c_diag = matrix.fullDiag();
            if (!matrix.isPrepared())
            {
                // copy upper band
                c_star[0] = c_star[0] / matrix.diag(0);
                for (int i = 1; i < size - 1; ++i)
                {
                    c_star[i] = c_star[i] / (matrix.diag(i) - c_star[i - 1] * matrix.lower(i - 1));
                }

                for (int i = 1; i < size; ++i)
                {
                    c_diag[i] = c_diag[i] - c_star[i - 1] * matrix.lower(i - 1);
                }
                matrix.setPrepared(true);
            }

            // MAIN LOOPS AFTER PREPARATION
            rhs[0] = rhs[0] / matrix.diag(0);
            for (int i = 1; i < size; ++i)
            {
                rhs[i] = (rhs[i] - rhs[i - 1] * matrix.lower(i - 1)) / c_diag[i];
            }

            // backtracking
            for (int i = size - 2; i >= 0; --i)
            {
                rhs[i] = rhs[i] - c_star[i] * rhs[i + 1];
            }
        }

        // does modify rhs but not the matrix
        inline void gaussPivotIP(const TDMatrix& matrix, std::vector<double>& rhs)
        {
            const auto n = matrix.size();
            assert(n >= 2);
            assert(n == rhs.size());

            // copy all bands
            auto du = matrix.fullUpper();
            auto d  = matrix.fullDiag();
            auto dl = matrix.fullLower();

            int i;
            for (i = 0; i < n - 2; ++i)
            {
                if (std::abs(d[i]) >= std::abs(dl[i]))
                {
                    // No row interchange required
                    // Caution: We assume here d[...] != 0
                    double fact = dl[i] / d[i];
                    d[i + 1] -= fact * du[i];
                    rhs[i + 1] -= fact * rhs[i];
                    dl[i] = 0.;
                }
                else
                {
                    // Interchange rows i and i+1
                    double fact = d[i] / dl[i];
                    d[i]        = dl[i];
                    double temp = d[i + 1];
                    d[i + 1]    = du[i] - fact * temp;
                    dl[i]       = du[i + 1];
                    du[i + 1]   = -fact * dl[i];
                    du[i]       = temp;
                    temp        = rhs[i];
                    rhs[i]      = rhs[i + 1];
                    rhs[i + 1]  = temp - fact * rhs[i + 1];
                }
            }

            // We have i = n - 2;
            if (std::abs(d[i]) >= std::abs(dl[i]))
            {
                // Caution: We assume here d[...] != 0
                double fact = dl[i] / d[i];
                d[i + 1] -= fact * du[i];
                rhs[i + 1] -= fact * rhs[i];
            }
            else
            {
                double fact = d[i] / dl[i];
                d[i]        = dl[i];
                double temp = d[i + 1];
                d[i + 1]    = du[i] - fact * temp;
                du[i]       = temp;
                temp        = rhs[i];
                rhs[i]      = rhs[i + 1];
                rhs[i + 1]  = temp - fact * rhs[i + 1];
            }

            // Back solve with the matrix U from the factorization.
            rhs[n - 1] /= d[n - 1];
            rhs[n - 2] = (rhs[n - 2] - du[n - 2] * rhs[n - 1]) / d[n - 2];
            for (i = n - 3; i >= 0; --i)
            {
                rhs[i] = (rhs[i] - du[i] * rhs[i + 1] - dl[i] * rhs[i + 2]) / d[i];
            }
        }

        // Modifies matrix (prepares matrix if necessary) and rhs
        inline void gaussReUsePivotIP(TDMatrix& matrix, std::vector<double>& rhs)
        {
            const auto n = matrix.size();
            assert(n >= 2);
            assert(n == rhs.size());

            // No copy, just references for convenience
            std::vector<double>& du  = matrix.fullUpper();
            std::vector<double>& d   = matrix.fullDiag();
            std::vector<double>& dl  = matrix.fullLower();
            std::vector<double>& du2 = matrix.fullSuperUpper();
            std::vector<int>& ipiv   = matrix.fullPivotIndex();

            if (!matrix.isPrepared())
            {
                // We have to prepare the matrix
                du2.resize(n - 2);
                ipiv.resize(n);

                // Initialize ipiv[i] = i and du2[i] = 0
                for (int i = 0; i < n; ++i)
                {
                    ipiv[i] = i;
                }
                for (int i = 0; i < n - 2; ++i)
                {
                    du2[i] = 0.;
                }
                int i = 0;
                for (; i < n - 2; ++i)
                {
                    if (std::abs(d[i]) >= std::abs(dl[i]))
                    {
                        // No row interchange required, eliminate dl[i]
                        // Caution: We assume d[i] != 0.
                        double fact = dl[i] / d[i];
                        dl[i]       = fact;
                        d[i + 1] -= fact * du[i];
                    }
                    else
                    {
                        // Interchange rows i and i+1, eliminate dl[i]
                        double fact = d[i] / dl[i];
                        d[i]        = dl[i];
                        dl[i]       = fact;
                        double temp = du[i];
                        du[i]       = d[i + 1];
                        d[i + 1]    = temp - fact * d[i + 1];
                        du2[i]      = du[i + 1];
                        du[i + 1]   = -fact * du[i + 1];
                        ipiv[i]     = i + 1;
                    }
                }
                // We have i = n - 2;
                if (std::abs(d[i]) >= std::abs(dl[i]))
                {
                    // Caution: We assume d[i] != 0.
                    double fact = dl[i] / d[i];
                    dl[i]       = fact;
                    d[i + 1] -= fact * du[i];
                }
                else
                {
                    double fact = d[i] / dl[i];
                    d[i]        = dl[i];
                    dl[i]       = fact;
                    double temp = du[i];
                    du[i]       = d[i + 1];
                    d[i + 1]    = temp - fact * d[i + 1];
                    ipiv[i]     = i + 1;
                }

                // Note: we assume that all d[...] != 0.
                matrix.setPrepared(true);
            }

            // Solve A*X = B using the LU factorization of A,
            // overwriting each right hand side vector with its solution.

            // Solve L*x = b.

            for (int i = 0; i < n - 1; ++i)
            {
                int ip      = ipiv[i];
                double temp = rhs[2 * i + 1 - ip] - dl[i] * rhs[ip];
                rhs[i]      = rhs[ip];
                rhs[i + 1]  = temp;
            }

            // Solve U*x = b.

            // We assume d[...] != 0.
            rhs[n - 1] /= d[n - 1];
            rhs[n - 2] = (rhs[n - 2] - du[n - 2] * rhs[n - 1]) / d[n - 2];
            for (int i = n - 3; i >= 0; --i)
            {
                rhs[i] = (rhs[i] - du[i] * rhs[i + 1] - du2[i] * rhs[i + 2]) / d[i];
            }
        }
    }
}

#endif  // ALGORITHMS_H
