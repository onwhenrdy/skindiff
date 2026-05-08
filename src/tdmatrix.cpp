#include "tdmatrix.h"

#include <algorithm>
#include <cmath>

namespace sc
{
    double TDMatrix::absMax() const noexcept
    {
        if (m_size < 1) return 0.0;

        double max_el = 0.0;
        for (int i = 0; i < m_size - 1; ++i)
        {
            max_el = std::max({max_el, std::abs(m_diag[i]), std::abs(m_lower[i]),
                               std::abs(m_upper[i])});
        }
        max_el = std::max(max_el, std::abs(m_diag[m_size - 1]));
        return max_el;
    }

    bool TDMatrix::isDiagonalDominant() const noexcept
    {
        for (int i = 1; i < m_size - 1; ++i)
        {
            if (m_diag[i] < (m_upper[i] + m_lower[i - 1]))
            {
                return false;
            }
        }
        return true;
    }

    void TDMatrix::multiplyBy(double val) noexcept
    {
        for (auto& d : m_diag)  d *= val;
        for (auto& d : m_upper) d *= val;
        for (auto& d : m_lower) d *= val;
    }
}
