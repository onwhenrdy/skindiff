#ifndef SC_MATRIXBUILDER_H
#define SC_MATRIXBUILDER_H

#include "compartment.h"
#include "geometry.h"
#include "sink.h"
#include "tdmatrix.h"
#include <functional>

namespace sc
{
    class MatrixBuilder
    {
      public:
        enum class Method
        {
            DSkin_1_3,  // central element cocentrations with back-flux correction
            DSkin_1_4,  // element edge concentrations (Crank MoD)
            DSkin_1_5   // fast version of 1_4
        };

        explicit MatrixBuilder(Method method = Method::DSkin_1_3);
        bool buildMatrix(const std::vector<Compartment>& compartments, const Geometry& geometry,
                         Sink* sink = nullptr);

        Method method() const;
        void setMethod(const Method& method);

        double maxModule() const;
        void setMaxModule(double max_module);

        const TDMatrix& matrixRhs() const;
        const TDMatrix& matrixLhs() const;
        int timesteps() const;

      private:
        bool build_M_DS_1_3(const std::vector<Compartment>& compartments, const Geometry& geometry,
                            Sink* sink = nullptr);

        bool build_M_DS_1_4(const std::vector<Compartment>& compartments, const Geometry& geometry,
                            Sink* sink = nullptr);

        bool build_M_DS_1_5(const std::vector<Compartment>& compartments, const Geometry& geometry,
                            Sink* sink = nullptr);

        // helper
        std::vector<double> createParamVector(int size,
                                              const std::vector<Compartment>& compartments,
                                              std::function<double(const Compartment&)> fun,
                                              Sink* sink = nullptr);

        // build the arithm. mean of values vec[i] and vec[j]
        double avgFromIdx(const std::vector<double>& vec, int i, int j);
        // build the harmonic mean of values vec[i] and vec[j]
        double harmMeanFromIdx(const std::vector<double>& vec, int i, int j);

        // calculate the back flux correction parameters k1-k4 from partition vector K for position idx
        void backFluxCorrection(const std::vector<double>& K, int idx, double& k1, double& k2,
                                double& k3, double& k4);
        // calculate the area correction parameters v1 and v2 from the area vector A for position idx
        void areaCorrection(const std::vector<double>& A, int idx, double& v1, double& v2);

        // Crank-Nicolson matrix for the LHS produced by matrix from the RHS of the equation
        TDMatrix fromRhs(const TDMatrix& matrix);

      private:
        Method m_method;
        double m_max_module;
        TDMatrix m_matrix_rhs;
        TDMatrix m_matrix_lhs;
        int m_timesteps;
    };

    inline double MatrixBuilder::avgFromIdx(const std::vector<double>& vec, int i, int j)
    {
        return 0.5 * (vec[i] + vec[j]);
    }

    inline double MatrixBuilder::harmMeanFromIdx(const std::vector<double>& vec, int i, int j)
    {
        if (vec[i] == vec[j])
        {
            return vec[i];
        }

        return 2.0 * vec[i] * vec[j] / (vec[i] + vec[j]);
    }

    inline void MatrixBuilder::backFluxCorrection(const std::vector<double>& K, int idx, double& k1,
                                                  double& k2, double& k3, double& k4)
    {
        if (K[idx + 1] > K[idx])
        {
            k2 = K[idx] / K[idx + 1];
        }
        else
        {
            k4 = K[idx + 1] / K[idx];
        }

        if (K[idx - 1] > K[idx])
        {
            k1 = K[idx] / K[idx - 1];
        }
        else
        {
            k3 = K[idx - 1] / K[idx];
        }
    }

    inline void MatrixBuilder::areaCorrection(const std::vector<double>& A, int idx, double& v1,
                                              double& v2)
    {
        // check if right compartment is smaller than xi
        if (A[idx + 1] < A[idx])
        {
            v1 = A[idx + 1] / A[idx];
        }

        // check if left compartment is smaller than xi
        if (A[idx - 1] < A[idx])
        {
            v2 = A[idx - 1] / A[idx];
        }
    }

    std::string toString(MatrixBuilder::Method method);
    MatrixBuilder::Method mbMethodFromString(const std::string& str, bool* ok = nullptr);
}
#endif  // MATRIXBUILDER_H
