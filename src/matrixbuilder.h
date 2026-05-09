#ifndef SC_MATRIXBUILDER_H
#define SC_MATRIXBUILDER_H

#include "compartment.h"
#include "geometry.h"
#include "sink.h"
#include "tdmatrix.h"

namespace sc
{
    // Cell-centred finite-volume builder in the activity variable u = c/K.
    // u is continuous at partition interfaces, so harmonic-mean face
    // conductances handle them naturally. Symmetric tri-diagonal (modulo
    // the absorbing sink BC), second-order in the interior. Time stepping
    // is Crank-Nicolson; the resulting LHS / RHS tri-diagonal matrices are
    // returned for the caller to use with the Thomas-reuse solver.
    class MatrixBuilder
    {
      public:
        MatrixBuilder() = default;

        bool buildMatrix(const std::vector<Compartment>& compartments,
                         const Geometry& geometry, Sink* sink = nullptr);

        [[nodiscard]] double maxModule() const noexcept { return m_max_module; }
        void setMaxModule(double max_module) noexcept { m_max_module = max_module; }

        [[nodiscard]] const TDMatrix& matrixRhs() const noexcept { return m_matrix_rhs; }
        [[nodiscard]] const TDMatrix& matrixLhs() const noexcept { return m_matrix_lhs; }
        [[nodiscard]] int timesteps() const noexcept { return m_timesteps; }

      private:
        double   m_max_module = 50.0;
        TDMatrix m_matrix_rhs;
        TDMatrix m_matrix_lhs;
        int      m_timesteps  = 1;
    };
}

#endif  // SC_MATRIXBUILDER_H
