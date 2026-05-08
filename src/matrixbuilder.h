#ifndef SC_MATRIXBUILDER_H
#define SC_MATRIXBUILDER_H

#include "compartment.h"
#include "geometry.h"
#include "sink.h"
#include "tdmatrix.h"

#include <optional>
#include <string>
#include <string_view>

namespace sc
{
    class MatrixBuilder
    {
      public:
        enum class Method
        {
            // Crank-style finite-difference in the c (concentration) variable.
            // Handles partition jumps via K_l/K_c factors in the stencil.
            // Roughly first-order at K interfaces.
            DSkin_1_4,

            // Cell-centred finite volume in the activity variable u = c/K.
            // u is continuous at partition interfaces, so harmonic-mean
            // face conductances handle them naturally. Symmetric tri-diagonal
            // (modulo the absorbing sink BC), second-order in the interior.
            Activity_FVM
        };

        // Whether `method` operates on activity (u = c/K) instead of c. Used
        // by System to convert the state vector at the simulation boundaries.
        [[nodiscard]] static bool isActivitySpace(Method m) noexcept
        {
            return m == Method::Activity_FVM;
        }

        explicit MatrixBuilder(Method method = Method::DSkin_1_4) noexcept : m_method(method) {}

        bool buildMatrix(const std::vector<Compartment>& compartments, const Geometry& geometry,
                         Sink* sink = nullptr);

        [[nodiscard]] Method method() const noexcept { return m_method; }
        void setMethod(Method method) noexcept { m_method = method; }

        [[nodiscard]] double maxModule() const noexcept { return m_max_module; }
        void setMaxModule(double max_module) noexcept { m_max_module = max_module; }

        [[nodiscard]] const TDMatrix& matrixRhs() const noexcept { return m_matrix_rhs; }
        [[nodiscard]] const TDMatrix& matrixLhs() const noexcept { return m_matrix_lhs; }
        [[nodiscard]] int timesteps() const noexcept { return m_timesteps; }

      private:
        Method   m_method     = Method::DSkin_1_4;
        double   m_max_module = 50.0;
        TDMatrix m_matrix_rhs;
        TDMatrix m_matrix_lhs;
        int      m_timesteps  = 1;
    };

    [[nodiscard]] std::string_view toString(MatrixBuilder::Method method) noexcept;
    [[nodiscard]] std::optional<MatrixBuilder::Method>
        mbMethodFromString(std::string_view str) noexcept;
}

#endif  // SC_MATRIXBUILDER_H
