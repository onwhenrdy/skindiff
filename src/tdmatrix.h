#ifndef SC_TDMATRIX_H
#define SC_TDMATRIX_H

#include <ostream>
#include <vector>
#include <cassert>

namespace sc
{
    // Triagional band matrix class
    class TDMatrix
    {
      public:
        TDMatrix();
        explicit TDMatrix(int size);
        void clear();

        // from 0 .. size-1
        double& diag(int idx);
        double diag(int idx) const;
        std::vector<double>& fullDiag();
        const std::vector<double>& fullDiag() const;
        bool isDiagonalDominat() const;

        // from 0 .. size- 2
        double& lower(int idx);
        double lower(int idx) const;
        std::vector<double>& fullLower();
        const std::vector<double>& fullLower() const;

        // from 0 ... size - 2
        double& upper(int idx);
        double upper(int idx) const;
        std::vector<double>& fullUpper();
        const std::vector<double>& fullUpper() const;

        // from 0 ... size - 3
        double& superUpper(int idx);
        double superUpper(int idx) const;
        std::vector<double>& fullSuperUpper();
        const std::vector<double>& fullSuperUpper() const;

        // from 0 ... size - 1
        int& pivotIndex(int idx);
        int pivotIndex(int idx) const;
        std::vector<int>& fullPivotIndex();
        const std::vector<int>& fullPivotIndex() const;

        double max() const;
        double absMax() const;

        void multiplyBy(double val);

        std::vector<double> operator*(const std::vector<double>& vec) const;
        void inlineMultiply(std::vector<double>& vec) const;

        // see https://de.wikipedia.org/wiki/Tridiagonalmatrix for index structure starting
        // with (0,0)
        // i = row, j = column
        // Beware: accessing elements outside the matrix band has undefined behavior
        double& operator()(int i, int j);
        double operator()(int i, int j) const;
        friend std::ostream& operator<<(std::ostream& out, const TDMatrix& matrix);

        // size is defined as the size of the diagonal
        // => Matrix is size x size
        int size() const;

        // returns whether the matrix is in factorized LU form
        // to be used for efficiently solving the corresponding
        // system of linear equations
        bool isPrepared() const;
        void setPrepared(bool prep);

      private:
        std::vector<double> m_upper;
        std::vector<double> m_diag;
        std::vector<double> m_lower;
        std::vector<double> m_superupper;
        std::vector<int> m_ipivot;

        int m_size;
        bool m_prepared;
    };

    std::ostream& operator<<(std::ostream& out, const TDMatrix& matrix);

    inline std::vector<double> TDMatrix::operator*(const std::vector<double>& vec) const
    {
        assert(m_size == vec.size());
        assert(m_size > 1);
        std::vector<double> result(m_size);

        result[0] = vec[0] * m_diag[0] + vec[1] * m_upper[0];
        for (int i = 1; i < m_size - 1; ++i)
        {
            result[i] = m_lower[i - 1] * vec[i - 1] + m_diag[i] * vec[i] + m_upper[i] * vec[i + 1];
        }
        const auto idx = m_size - 1;
        result[idx]    = m_lower[idx - 1] * vec[idx - 1] + m_diag[idx] * vec[idx];

        return result;
    }

    inline void TDMatrix::inlineMultiply(std::vector<double> &vec) const
    {
        assert(m_size == vec.size());
        assert(m_size > 1);

        auto tmp = vec[0];
        vec[0] = vec[0] * m_diag[0] + vec[1] * m_upper[0];
        for (int i = 1; i < m_size - 1; ++i)
        {
            const auto old_val_i = vec[i];
            vec[i] = m_lower[i - 1] * tmp + m_diag[i] * old_val_i + m_upper[i] * vec[i + 1];
            tmp = old_val_i;
        }
        const auto idx = m_size - 1;
        vec[idx]    = m_lower[idx - 1] * tmp + m_diag[idx] * vec[idx];
    }
}

#endif  // TDMATRIX_H
