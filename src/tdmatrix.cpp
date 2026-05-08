#include "tdmatrix.h"
#include <algorithm>
#include <ostream>
#include <sstream>

namespace sc
{
    TDMatrix::TDMatrix() : m_size(0), m_prepared(false)
    {
    }

    TDMatrix::TDMatrix(int size) : m_size(size), m_prepared(false)
    {
        m_diag.resize(m_size, 0.0);
        m_lower.resize(m_size - 1, 0.0);
        m_upper.resize(m_size - 1, 0.0);
    }

    void TDMatrix::clear()
    {
        for (auto& val : m_diag)
        {
            val = 0.0;
        }

        for (auto& val : m_upper)
        {
            val = 0.0;
        }

        for (auto& val : m_lower)
        {
            val = 0.0;
        }
    }

    double& TDMatrix::diag(int idx)
    {
        assert(idx >= 0 && idx < m_size);
        return m_diag[idx];
    }

    double TDMatrix::diag(int idx) const
    {
        assert(idx >= 0 && idx < m_size);
        return m_diag[idx];
    }

    std::vector<double>& TDMatrix::fullDiag()
    {
        return m_diag;
    }

    const std::vector<double>& TDMatrix::fullDiag() const
    {
        return m_diag;
    }

    bool TDMatrix::isDiagonalDominat() const
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

    double& TDMatrix::lower(int idx)
    {
        assert(idx >= 0 && idx < m_size - 1);
        return m_lower[idx];
    }

    double TDMatrix::lower(int idx) const
    {
        assert(idx >= 0 && idx < m_size - 1);
        return m_lower[idx];
    }

    std::vector<double>& TDMatrix::fullLower()
    {
        return m_lower;
    }

    const std::vector<double>& TDMatrix::fullLower() const
    {
        return m_lower;
    }

    double& TDMatrix::upper(int idx)
    {
        assert(idx >= 0 && idx < m_size - 1);
        return m_upper[idx];
    }

    double TDMatrix::upper(int idx) const
    {
        assert(idx >= 0 && idx < m_size - 1);
        return m_upper[idx];
    }

    std::vector<double>& TDMatrix::fullUpper()
    {
        return m_upper;
    }

    const std::vector<double>& TDMatrix::fullUpper() const
    {
        return m_upper;
    }

    double& TDMatrix::superUpper(int idx)
    {
        assert(idx >= 0 && idx < m_size - 2);
        return m_superupper[idx];
    }

    double TDMatrix::superUpper(int idx) const
    {
        assert(idx >= 0 && idx < m_size - 2);
        return m_superupper[idx];
    }

    std::vector<double>& TDMatrix::fullSuperUpper()
    {
        return m_superupper;
    }

    const std::vector<double>& TDMatrix::fullSuperUpper() const
    {
        return m_superupper;
    }

    int& TDMatrix::pivotIndex(int idx)
    {
        assert(idx >= 0 && idx < m_size - 2);
        return m_ipivot[idx];
    }

    int TDMatrix::pivotIndex(int idx) const
    {
        assert(idx >= 0 && idx < m_size - 2);
        return m_ipivot[idx];
    }

    std::vector<int>& TDMatrix::fullPivotIndex()
    {
        return m_ipivot;
    }

    const std::vector<int>& TDMatrix::fullPivotIndex() const
    {
        return m_ipivot;
    }

    double TDMatrix::max() const
    {
        auto max_el = 0.0;
        if (m_size < 1)
        {
            return max_el;
        }

        for (int i = 0; i < m_size - 1; ++i)
        {
            max_el = std::max({max_el, m_diag[i], m_lower[i], m_upper[i]});
        }
        max_el = std::max(max_el, m_diag[m_size - 1]);

        return max_el;
    }

    double TDMatrix::absMax() const
    {
        auto max_el = 0.0;
        if (m_size < 1)
        {
            return max_el;
        }

        for (int i = 0; i < m_size - 1; ++i)
        {
            max_el = std::max({std::abs(max_el), std::abs(m_diag[i]), std::abs(m_lower[i]),
                               std::abs(m_upper[i])});
        }
        max_el = std::max(max_el, std::abs(m_diag[m_size - 1]));

        return max_el;
    }

    void TDMatrix::multiplyBy(double val)
    {
        for (auto& d : m_diag)
        {
            d *= val;
        }

        for (auto& d : m_upper)
        {
            d *= val;
        }

        for (auto& d : m_lower)
        {
            d *= val;
        }
    }

    // row, col
    double& TDMatrix::operator()(int i, int j)
    {
        assert(i < m_size && j < m_size);
        assert(i >= 0 && j >= 0);
        assert(j <= i + 1 || j >= i - 1);

        if (i == j)
        {
            return m_diag[i];
        }
        // row < column -> upper
        else if (i < j)
        {
            return m_upper[i];
        }

        return m_lower[i - 1];
    }

    double TDMatrix::operator()(int i, int j) const
    {
        assert(i < m_size && j < m_size);
        assert(i >= 0 && j >= 0);
        assert(j <= i + 1 || j >= i - 1);

        if (i == j)
        {
            return m_diag[i];
        }
        // row < column -> upper
        else if (i < j)
        {
            return m_upper[i];
        }

        return m_lower[i - 1];
    }

    int TDMatrix::size() const
    {
        return m_size;
    }

    bool TDMatrix::isPrepared() const
    {
        return m_prepared;
    }

    void TDMatrix::setPrepared(bool prep)
    {
        m_prepared = prep;
    }

    std::ostream& operator<<(std::ostream& out, const TDMatrix& matrix)
    {
        auto vec_out = [](const std::vector<double>& vec) {
            std::stringstream ss;
            for (auto val : vec)
            {
                ss << val << " ";
            }
            return ss.str();
        };

        out << "Size : " << matrix.size() << " diagonal elements\n";
        out << "Prep : " << (matrix.isPrepared() ? "yes" : "no") << "\n";
        out << "UDEs : " << vec_out(matrix.fullUpper()) << "\n";
        out << "CDEs : " << vec_out(matrix.fullLower()) << "\n";
        out << "LDEs : " << vec_out(matrix.fullLower()) << "\n";

        return out;
    }
}
