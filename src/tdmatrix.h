#ifndef SC_TDMATRIX_H
#define SC_TDMATRIX_H

#include <cassert>
#include <cstddef>
#include <vector>

namespace sc
{
    // Tri-diagonal band matrix.
    //
    // Indexing follows row-major notation, see
    //     https://de.wikipedia.org/wiki/Tridiagonalmatrix
    // for the index structure starting with (0, 0).
    //
    // After a Thomas-style solve, the matrix may be in factored form;
    // isPrepared() / setPrepared() track that state.
    class TDMatrix
    {
      public:
        TDMatrix() = default;
        explicit TDMatrix(int size)
            : m_diag(size, 0.0), m_lower(size > 0 ? size - 1 : 0, 0.0),
              m_upper(size > 0 ? size - 1 : 0, 0.0), m_size(size)
        {
        }

        // size of the diagonal (matrix is size x size)
        [[nodiscard]] int size() const noexcept { return m_size; }

        double& diag(int i)
        {
            assert(i >= 0 && i < m_size);
            return m_diag[i];
        }
        double diag(int i) const
        {
            assert(i >= 0 && i < m_size);
            return m_diag[i];
        }
        std::vector<double>& fullDiag() noexcept { return m_diag; }
        const std::vector<double>& fullDiag() const noexcept { return m_diag; }

        double& lower(int i)
        {
            assert(i >= 0 && i < m_size - 1);
            return m_lower[i];
        }
        double lower(int i) const
        {
            assert(i >= 0 && i < m_size - 1);
            return m_lower[i];
        }
        std::vector<double>& fullLower() noexcept { return m_lower; }
        const std::vector<double>& fullLower() const noexcept { return m_lower; }

        double& upper(int i)
        {
            assert(i >= 0 && i < m_size - 1);
            return m_upper[i];
        }
        double upper(int i) const
        {
            assert(i >= 0 && i < m_size - 1);
            return m_upper[i];
        }
        std::vector<double>& fullUpper() noexcept { return m_upper; }
        const std::vector<double>& fullUpper() const noexcept { return m_upper; }

        [[nodiscard]] double absMax() const noexcept;
        [[nodiscard]] bool isDiagonalDominant() const noexcept;

        void multiplyBy(double val) noexcept;

        // y = M * x
        [[nodiscard]] std::vector<double> operator*(const std::vector<double>& vec) const;
        // x = M * x (in place)
        void inlineMultiply(std::vector<double>& vec) const;

        [[nodiscard]] bool isPrepared() const noexcept { return m_prepared; }
        void setPrepared(bool prep) noexcept { m_prepared = prep; }

      private:
        std::vector<double> m_diag;
        std::vector<double> m_lower;
        std::vector<double> m_upper;
        int  m_size     = 0;
        bool m_prepared = false;
    };

    inline std::vector<double> TDMatrix::operator*(const std::vector<double>& vec) const
    {
        assert(static_cast<std::size_t>(m_size) == vec.size());
        assert(m_size > 1);

        std::vector<double> result(static_cast<std::size_t>(m_size));
        result[0] = vec[0] * m_diag[0] + vec[1] * m_upper[0];
        for (int i = 1; i < m_size - 1; ++i)
        {
            result[i] = m_lower[i - 1] * vec[i - 1] + m_diag[i] * vec[i] + m_upper[i] * vec[i + 1];
        }
        const auto idx = m_size - 1;
        result[idx]    = m_lower[idx - 1] * vec[idx - 1] + m_diag[idx] * vec[idx];
        return result;
    }

    inline void TDMatrix::inlineMultiply(std::vector<double>& vec) const
    {
        assert(static_cast<std::size_t>(m_size) == vec.size());
        assert(m_size > 1);

        auto tmp = vec[0];
        vec[0]   = vec[0] * m_diag[0] + vec[1] * m_upper[0];
        for (int i = 1; i < m_size - 1; ++i)
        {
            const auto old_val_i = vec[i];
            vec[i] = m_lower[i - 1] * tmp + m_diag[i] * old_val_i + m_upper[i] * vec[i + 1];
            tmp    = old_val_i;
        }
        const auto idx = m_size - 1;
        vec[idx]       = m_lower[idx - 1] * tmp + m_diag[idx] * vec[idx];
    }
}

#endif  // SC_TDMATRIX_H
