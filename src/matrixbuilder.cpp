#include "matrixbuilder.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <vector>

namespace sc
{
    namespace
    {
        // Build a per-cell parameter vector by broadcasting the compartment-level
        // value into every cell that belongs to it. Cells covered by the sink
        // inherit the value of the cell immediately above (the last skin cell).
        template <typename Fun>
        std::vector<double>
        createParamVector(int size, const std::vector<Compartment>& compartments, Fun fun,
                          const Sink* sink = nullptr)
        {
            std::vector<double> result(static_cast<std::size_t>(size), 0.0);

            for (const auto& comp : compartments)
            {
                const auto val = fun(comp);
                for (int i = comp.geo_from; i <= comp.geo_to; ++i)
                {
                    result[static_cast<std::size_t>(i)] = val;
                }
            }

            if (sink)
            {
                const auto idx = sink->geo_from;
                if (idx > 0)
                {
                    result[static_cast<std::size_t>(idx)] =
                        result[static_cast<std::size_t>(idx - 1)];
                }
            }

            return result;
        }
    }

    bool MatrixBuilder::buildMatrix(const std::vector<Compartment>& compartments,
                                    const Geometry& geometry, Sink* sink)
    {
        assert(m_method == Method::DSkin_1_4);
        assert(!compartments.empty());

        const auto sys_size = geometry.size();
        const auto& space_steps = geometry.spaceSteps();

        const auto D_vec = createParamVector(sys_size, compartments,
                                             [](const Compartment& c) { return c.D; }, sink);
        const auto K_vec = createParamVector(sys_size, compartments,
                                             [](const Compartment& c) { return c.K; }, sink);
        const auto A_vec = createParamVector(sys_size, compartments,
                                             [](const Compartment& c) { return c.area_um2; }, sink);

        m_matrix_rhs = TDMatrix(sys_size);
        m_matrix_lhs = TDMatrix(sys_size);

        // Reflecting boundary at the start.
        {
            const auto l_c = space_steps[0];
            const auto l_r = space_steps[1];
            const auto D_c = D_vec[0];
            const auto D_r = D_vec[1];
            const auto K_c = K_vec[0];
            const auto K_r = K_vec[1];
            const auto A_c = A_vec[0];
            const auto A_r = A_vec[1];

            const auto h2 = (l_c + l_r) / 2.0;
            const auto upper_f =
                (l_c + l_r) * D_c * D_r / (l_c * D_r + K_c / K_r * l_r * D_c) / (h2 * h2);

            const auto upper_val = upper_f * K_c / K_r * std::min(1.0, A_r / A_c);
            const auto mid_val   = upper_f * std::min(1.0, A_r / A_c);

            m_matrix_rhs.diag(0)  = mid_val;
            m_matrix_rhs.upper(0) = upper_val;
        }

        // Main building loop.
        for (int i = 1; i < sys_size - 1; ++i)
        {
            const auto l_l = space_steps[static_cast<std::size_t>(i - 1)];
            const auto l_c = space_steps[static_cast<std::size_t>(i)];
            const auto l_r = space_steps[static_cast<std::size_t>(i + 1)];

            const auto D_l = D_vec[static_cast<std::size_t>(i - 1)];
            const auto D_c = D_vec[static_cast<std::size_t>(i)];
            const auto D_r = D_vec[static_cast<std::size_t>(i + 1)];

            const auto K_l = K_vec[static_cast<std::size_t>(i - 1)];
            const auto K_c = K_vec[static_cast<std::size_t>(i)];
            const auto K_r = K_vec[static_cast<std::size_t>(i + 1)];

            const auto A_l = A_vec[static_cast<std::size_t>(i - 1)];
            const auto A_c = A_vec[static_cast<std::size_t>(i)];
            const auto A_r = A_vec[static_cast<std::size_t>(i + 1)];

            const auto h1 = (l_l + l_c) / 2.0;
            const auto h2 = (l_c + l_r) / 2.0;

            // Weighted harmonic-mean coefficients corrected for partition and area.
            const auto lower_f = (l_l + l_c) * D_l * D_c / (l_l * D_c + K_l / K_c * l_c * D_l) *
                                 2.0 * h2 / (h1 * h2 * (h1 + h2));
            const auto upper_f = (l_c + l_r) * D_c * D_r / (l_c * D_r + K_c / K_r * l_r * D_c) *
                                 2.0 * h1 / (h1 * h2 * (h1 + h2));

            const auto lower_val = lower_f * std::min(1.0, A_l / A_c);
            const auto upper_val = upper_f * K_c / K_r * std::min(1.0, A_r / A_c);
            const auto mid_val   = lower_f * K_l / K_c * std::min(1.0, A_l / A_c) +
                                   upper_f * std::min(1.0, A_r / A_c);

            m_matrix_rhs.diag(i)      = mid_val;
            m_matrix_rhs.upper(i)     = upper_val;
            m_matrix_rhs.lower(i - 1) = lower_val;
        }

        // Bottom row, sink-side reflecting term.
        {
            const auto l_l = space_steps[static_cast<std::size_t>(sys_size - 2)];
            const auto l_c = space_steps[static_cast<std::size_t>(sys_size - 1)];
            const auto D_l = D_vec[static_cast<std::size_t>(sys_size - 2)];
            const auto D_c = D_vec[static_cast<std::size_t>(sys_size - 1)];
            const auto K_l = K_vec[static_cast<std::size_t>(sys_size - 2)];
            const auto K_c = K_vec[static_cast<std::size_t>(sys_size - 1)];
            const auto A_l = A_vec[static_cast<std::size_t>(sys_size - 2)];
            const auto A_c = A_vec[static_cast<std::size_t>(sys_size - 1)];

            const auto h1 = (l_l + l_c) / 2.0;
            const auto lower_f =
                (l_l + l_c) * D_l * D_c / (l_l * D_c + K_l / K_c * l_c * D_l) / (h1 * h1);

            const auto lower_val = lower_f * std::min(1.0, A_l / A_c);
            const auto mid_val   = lower_f * K_l / K_c * std::min(1.0, A_l / A_c);

            m_matrix_rhs.diag(sys_size - 1)  = mid_val;
            m_matrix_rhs.lower(sys_size - 2) = lower_val;
        }

        // Apply dt scaling and split into Crank-Nicolson LHS / RHS.
        const auto max_m = m_matrix_rhs.absMax();
        m_timesteps      = static_cast<int>(std::max(1.0, std::ceil(max_m / m_max_module)));
        const auto dt    = 1.0 / m_timesteps;
        m_matrix_rhs.multiplyBy(dt);

        for (int i = 0; i < sys_size - 1; ++i)
        {
            m_matrix_lhs.diag(i)  = 2.0 + m_matrix_rhs.diag(i);
            m_matrix_lhs.lower(i) = -m_matrix_rhs.lower(i);
            m_matrix_lhs.upper(i) = -m_matrix_rhs.upper(i);

            m_matrix_rhs.diag(i) = 2.0 - m_matrix_rhs.diag(i);
        }

        // No upward flux from the sink cell.
        m_matrix_rhs.upper(sys_size - 2) = 0.0;
        m_matrix_lhs.upper(sys_size - 2) = 0.0;

        // Sink behaviour.
        if (sink)
        {
            m_matrix_rhs.diag(sys_size - 1) = 2.0;
            m_matrix_lhs.diag(sys_size - 1) = 2.0;

            if (sink->type == Sink::Type::PK_Compartment)
            {
                m_matrix_rhs.diag(sys_size - 1) = 2.0 - dt * sink->kEl();
                m_matrix_lhs.diag(sys_size - 1) = 2.0 + dt * sink->kEl();
            }
        }

        // Infinite-dose donor (concentration clamped at the top boundary).
        if (!compartments.front().finite_dose)
        {
            m_matrix_rhs.diag(0) = 2.0;
            m_matrix_lhs.diag(0) = 2.0;

            m_matrix_rhs.upper(0) = 0.0;
            m_matrix_lhs.upper(0) = 0.0;
        }

        return true;
    }

    std::string_view toString(MatrixBuilder::Method method) noexcept
    {
        switch (method)
        {
            case MatrixBuilder::Method::DSkin_1_4: return "DSkin_1_4";
        }
        return "DSkin_1_4";
    }

    std::optional<MatrixBuilder::Method> mbMethodFromString(std::string_view str) noexcept
    {
        if (str == "DSkin_1_4") return MatrixBuilder::Method::DSkin_1_4;
        return std::nullopt;
    }
}
