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

    namespace
    {
        bool buildActivityFVM(double max_module, int& out_timesteps,
                              TDMatrix& out_rhs, TDMatrix& out_lhs,
                              const std::vector<Compartment>& compartments,
                              const Geometry& geometry, Sink* sink);
    }

    bool MatrixBuilder::buildMatrix(const std::vector<Compartment>& compartments,
                                    const Geometry& geometry, Sink* sink)
    {
        assert(!compartments.empty());

        if (m_method == Method::Activity_FVM)
        {
            return buildActivityFVM(m_max_module, m_timesteps, m_matrix_rhs, m_matrix_lhs,
                                    compartments, geometry, sink);
        }

        assert(m_method == Method::DSkin_1_4);

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

        // Infinite-dose donor: hold every donor cell at its initial value.
        // (The previous implementation only clamped cell 0, which left donor
        // cells 1..N free to drift slightly during transients. Clamping the
        // full compartment makes the donor act as a true Dirichlet reservoir
        // at c_init throughout, regardless of D_donor and donor mesh density.)
        if (!compartments.front().finite_dose)
        {
            const auto& donor = compartments.front();
            for (int i = donor.geo_from; i <= donor.geo_to; ++i)
            {
                m_matrix_rhs.diag(i) = 2.0;
                m_matrix_lhs.diag(i) = 2.0;
                if (i > 0)
                {
                    m_matrix_rhs.lower(i - 1) = 0.0;
                    m_matrix_lhs.lower(i - 1) = 0.0;
                }
                if (i < sys_size - 1)
                {
                    m_matrix_rhs.upper(i) = 0.0;
                    m_matrix_lhs.upper(i) = 0.0;
                }
            }
        }

        return true;
    }

    std::string_view toString(MatrixBuilder::Method method) noexcept
    {
        switch (method)
        {
            case MatrixBuilder::Method::DSkin_1_4:    return "DSkin_1_4";
            case MatrixBuilder::Method::Activity_FVM: return "Activity_FVM";
        }
        return "DSkin_1_4";
    }

    std::optional<MatrixBuilder::Method> mbMethodFromString(std::string_view str) noexcept
    {
        if (str == "DSkin_1_4")    return MatrixBuilder::Method::DSkin_1_4;
        if (str == "Activity_FVM") return MatrixBuilder::Method::Activity_FVM;
        return std::nullopt;
    }

    // ===========================================================================
    // Activity-FVM matrix builder (cell-centred finite volume in u = c/K).
    //
    // Equation: theta_i * h_i * du_i/dt = alpha_l (u_{i-1} - u_i) + alpha_r (u_{i+1} - u_i)
    // where:
    //   theta = K * A         (cell capacity weight)
    //   kappa = K * A * D     (face conductance weight)
    //   alpha_face = 2 * kappa_l * kappa_r / (h_l * kappa_r + h_r * kappa_l)
    //
    // The bulk operator is symmetric. The membrane <-> sink boundary is treated
    // by (a) using the cell-edge Dirichlet limit (kappa_sink -> infinity) so the
    // face conductance there reduces to 2 * kappa_membrane / h_membrane, and
    // (b) zeroing the membrane row's upper coupling so the membrane sees u=0 at
    // the sink while the sink cell continues to accumulate the throughput.
    // This matches the DSkin_1_4 absorbing-BC convention for back-compat.
    // ===========================================================================
    namespace
    {
        bool buildActivityFVM(double max_module, int& out_timesteps,
                              TDMatrix& out_rhs, TDMatrix& out_lhs,
                              const std::vector<Compartment>& compartments,
                              const Geometry& geometry, Sink* sink)
        {
            const auto N         = geometry.size();
            const auto& h        = geometry.spaceSteps();
            assert(N > 1);

            // Per-cell K, A, D from compartment broadcasting (sink inherits from
            // last skin layer, by createParamVector convention).
            auto K_vec = createParamVector(N, compartments,
                                           [](const Compartment& c) { return c.K; }, sink);
            auto D_vec = createParamVector(N, compartments,
                                           [](const Compartment& c) { return c.D; }, sink);
            auto A_vec = createParamVector(N, compartments,
                                           [](const Compartment& c) { return c.area_um2; }, sink);

            // Sink cell convention: K_sink = 1 (the activity in the receiver
            // compartment IS its concentration).
            if (sink)
            {
                K_vec[static_cast<std::size_t>(sink->geo_from)] = 1.0;
            }

            std::vector<double> kappa(static_cast<std::size_t>(N));
            std::vector<double> theta(static_cast<std::size_t>(N));
            for (int i = 0; i < N; ++i)
            {
                kappa[static_cast<std::size_t>(i)] =
                    K_vec[static_cast<std::size_t>(i)] *
                    A_vec[static_cast<std::size_t>(i)] *
                    D_vec[static_cast<std::size_t>(i)];
                theta[static_cast<std::size_t>(i)] =
                    K_vec[static_cast<std::size_t>(i)] *
                    A_vec[static_cast<std::size_t>(i)];
            }

            // Face conductances alpha_{i+1/2} for i = 0..N-2.
            std::vector<double> alpha(static_cast<std::size_t>(N - 1), 0.0);
            for (int i = 0; i < N - 1; ++i)
            {
                const auto k_l = kappa[static_cast<std::size_t>(i)];
                const auto k_r = kappa[static_cast<std::size_t>(i + 1)];
                const auto h_l = h[static_cast<std::size_t>(i)];
                const auto h_r = h[static_cast<std::size_t>(i + 1)];
                const auto den = h_l * k_r + h_r * k_l;
                alpha[static_cast<std::size_t>(i)] = (den > 0) ? 2.0 * k_l * k_r / den : 0.0;
            }

            // Membrane <-> sink face: treat sink as Dirichlet at the cell edge
            // (kappa_sink -> infinity in the harmonic mean), giving the half-cell
            // conductance of the membrane bottom: 2 * kappa_mem / h_mem.
            if (sink)
            {
                const auto mem_idx = sink->geo_from - 1;
                if (mem_idx >= 0)
                {
                    const auto k_mem = kappa[static_cast<std::size_t>(mem_idx)];
                    const auto h_mem = h[static_cast<std::size_t>(mem_idx)];
                    alpha[static_cast<std::size_t>(mem_idx)] = 2.0 * k_mem / h_mem;
                }
            }

            // Donor <-> first-skin-layer face for an infinite-dose donor:
            // treat the donor as Dirichlet at the cell edge (kappa_donor ->
            // infinity in the harmonic mean), giving the half-cell conductance
            // of the first skin layer cell: 2 * kappa_first_skin / h_first_skin.
            // This eliminates the residual "fast-but-finite donor" artefact in
            // the natural FVM coefficient.
            if (!compartments.front().finite_dose)
            {
                const auto& donor = compartments.front();
                if (donor.geo_to + 1 < N)
                {
                    const auto skin_idx = donor.geo_to + 1;
                    const auto k_skin = kappa[static_cast<std::size_t>(skin_idx)];
                    const auto h_skin = h[static_cast<std::size_t>(skin_idx)];
                    alpha[static_cast<std::size_t>(donor.geo_to)] = 2.0 * k_skin / h_skin;
                }
            }

            // Assemble |M| where M is the spatial operator
            //   theta_i * h_i * du_i/dt = -|M_diag| * u_i + |M_lower| * u_{i-1}
            //                                          + |M_upper| * u_{i+1}
            // Storage convention matches DSkin_1_4: positive magnitudes on every
            // band; CN sign-flip happens after dt scaling.
            out_rhs = TDMatrix(N);

            // Top boundary: reflecting (no left face).
            {
                const auto a_r   = alpha[0];
                const auto th_h0 = theta[0] * h[0];
                out_rhs.diag(0)  = a_r / th_h0;
                out_rhs.upper(0) = a_r / th_h0;
            }

            // Interior cells.
            for (int i = 1; i < N - 1; ++i)
            {
                const auto a_l   = alpha[static_cast<std::size_t>(i - 1)];
                const auto a_r   = alpha[static_cast<std::size_t>(i)];
                const auto th_hi = theta[static_cast<std::size_t>(i)] *
                                   h[static_cast<std::size_t>(i)];
                out_rhs.lower(i - 1) = a_l / th_hi;
                out_rhs.diag(i)      = (a_l + a_r) / th_hi;
                out_rhs.upper(i)     = a_r / th_hi;
            }

            // Last cell (sink in the simulator's setup): natural FVM closure with
            // no right-hand neighbour. We rewrite the diag below if a sink is
            // present.
            {
                const auto a_l   = alpha[static_cast<std::size_t>(N - 2)];
                const auto th_hi = theta[static_cast<std::size_t>(N - 1)] *
                                   h[static_cast<std::size_t>(N - 1)];
                out_rhs.lower(N - 2) = a_l / th_hi;
                out_rhs.diag(N - 1)  = a_l / th_hi;
            }

            // Pick dt / sub-step count from the largest |M| band entry.
            const auto max_m  = out_rhs.absMax();
            out_timesteps     = static_cast<int>(std::max(1.0, std::ceil(max_m / max_module)));
            const auto dt     = 1.0 / out_timesteps;
            out_rhs.multiplyBy(dt);

            // Crank-Nicolson sign-flip + symmetric LHS.
            out_lhs = TDMatrix(N);
            for (int i = 0; i < N - 1; ++i)
            {
                out_lhs.diag(i)  = 2.0 + out_rhs.diag(i);
                out_lhs.lower(i) = -out_rhs.lower(i);
                out_lhs.upper(i) = -out_rhs.upper(i);

                out_rhs.diag(i)  = 2.0 - out_rhs.diag(i);
            }

            // Sink BC: decouple the membrane row from the sink (no upward flux),
            // and reset the sink diag to the perfect / PK form.
            if (sink)
            {
                out_rhs.upper(N - 2) = 0.0;
                out_lhs.upper(N - 2) = 0.0;

                if (sink->type == Sink::Type::PK_Compartment)
                {
                    out_rhs.diag(N - 1) = 2.0 - dt * sink->kEl();
                    out_lhs.diag(N - 1) = 2.0 + dt * sink->kEl();
                }
                else
                {
                    out_rhs.diag(N - 1) = 2.0;
                    out_lhs.diag(N - 1) = 2.0;
                }
            }

            // Infinite-dose donor: hold every donor cell at its initial
            // value. Combined with the alpha override above, this makes the
            // donor act as a true Dirichlet reservoir at the donor-membrane
            // interface, regardless of D_donor or how many donor cells the
            // mesh produces.
            if (!compartments.front().finite_dose)
            {
                const auto& donor = compartments.front();
                for (int i = donor.geo_from; i <= donor.geo_to; ++i)
                {
                    out_rhs.diag(i) = 2.0;
                    out_lhs.diag(i) = 2.0;
                    if (i > 0)
                    {
                        out_rhs.lower(i - 1) = 0.0;
                        out_lhs.lower(i - 1) = 0.0;
                    }
                    if (i < N - 1)
                    {
                        out_rhs.upper(i) = 0.0;
                        out_lhs.upper(i) = 0.0;
                    }
                }
            }

            return true;
        }
    }
}
