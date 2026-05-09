#include "matrixbuilder.h"

#include <algorithm>
#include <cassert>
#include <cmath>
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

    // ===========================================================================
    // Cell-centred finite-volume builder in the activity variable u = c/K.
    //
    //   theta_i * h_i * du_i/dt = alpha_l (u_{i-1} - u_i) + alpha_r (u_{i+1} - u_i)
    //
    // with
    //   theta = K * A         (cell capacity weight)
    //   kappa = K * A * D     (face conductance weight)
    //   alpha_face = 2 * kappa_l * kappa_r / (h_l * kappa_r + h_r * kappa_l)
    //
    // The bulk operator is symmetric. Boundary handling:
    //
    //  - Top (donor surface): reflecting by default; for an infinite-dose
    //    donor we clamp every donor cell at its initial value AND set the
    //    donor / first-skin face conductance to its kappa_donor -> infinity
    //    limit (2 * kappa_skin / h_skin), giving a true Dirichlet BC.
    //
    //  - Bottom (membrane <-> sink): kappa_sink -> infinity in the harmonic
    //    mean (the sink is a Dirichlet phantom from the membrane's point of
    //    view). The membrane row's upper coupling to the sink is zeroed out;
    //    the sink cell still receives the throughput via its lower
    //    coefficient and serves as the cumulative-mass accumulator.
    // ===========================================================================
    bool MatrixBuilder::buildMatrix(const std::vector<Compartment>& compartments,
                                    const Geometry& geometry, Sink* sink)
    {
        assert(!compartments.empty());

        const auto N         = geometry.size();
        const auto& h        = geometry.spaceSteps();
        assert(N > 1);

        // Per-cell K, A, D from compartment broadcasting (sink inherits from
        // the last skin layer).
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

        // Membrane <-> sink face: cell-edge Dirichlet limit.
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

        // Donor <-> first-skin face for an infinite-dose donor: cell-edge
        // Dirichlet limit on the skin side.
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
        // Storage convention: positive magnitudes on every band; the CN
        // sign-flip is applied after the dt scaling below.
        m_matrix_rhs = TDMatrix(N);

        // Top boundary: reflecting (no left face).
        {
            const auto a_r   = alpha[0];
            const auto th_h0 = theta[0] * h[0];
            m_matrix_rhs.diag(0)  = a_r / th_h0;
            m_matrix_rhs.upper(0) = a_r / th_h0;
        }

        // Interior cells.
        for (int i = 1; i < N - 1; ++i)
        {
            const auto a_l   = alpha[static_cast<std::size_t>(i - 1)];
            const auto a_r   = alpha[static_cast<std::size_t>(i)];
            const auto th_hi = theta[static_cast<std::size_t>(i)] *
                               h[static_cast<std::size_t>(i)];
            m_matrix_rhs.lower(i - 1) = a_l / th_hi;
            m_matrix_rhs.diag(i)      = (a_l + a_r) / th_hi;
            m_matrix_rhs.upper(i)     = a_r / th_hi;
        }

        // Last cell (sink in the simulator's setup): natural FVM closure with
        // no right-hand neighbour. The diag is rewritten below if a sink is
        // present.
        {
            const auto a_l   = alpha[static_cast<std::size_t>(N - 2)];
            const auto th_hi = theta[static_cast<std::size_t>(N - 1)] *
                               h[static_cast<std::size_t>(N - 1)];
            m_matrix_rhs.lower(N - 2) = a_l / th_hi;
            m_matrix_rhs.diag(N - 1)  = a_l / th_hi;
        }

        // Pick dt / sub-step count from the largest |M| band entry.
        const auto max_m  = m_matrix_rhs.absMax();
        m_timesteps       = static_cast<int>(std::max(1.0, std::ceil(max_m / m_max_module)));
        const auto dt     = 1.0 / m_timesteps;
        m_matrix_rhs.multiplyBy(dt);

        // Crank-Nicolson sign-flip + symmetric LHS.
        m_matrix_lhs = TDMatrix(N);
        for (int i = 0; i < N - 1; ++i)
        {
            m_matrix_lhs.diag(i)  = 2.0 + m_matrix_rhs.diag(i);
            m_matrix_lhs.lower(i) = -m_matrix_rhs.lower(i);
            m_matrix_lhs.upper(i) = -m_matrix_rhs.upper(i);

            m_matrix_rhs.diag(i)  = 2.0 - m_matrix_rhs.diag(i);
        }

        // Sink BC: decouple the membrane row from the sink (no upward flux),
        // reset the sink diag to the perfect / PK form.
        if (sink)
        {
            m_matrix_rhs.upper(N - 2) = 0.0;
            m_matrix_lhs.upper(N - 2) = 0.0;

            if (sink->type == Sink::Type::PK_Compartment)
            {
                m_matrix_rhs.diag(N - 1) = 2.0 - dt * sink->kEl();
                m_matrix_lhs.diag(N - 1) = 2.0 + dt * sink->kEl();
            }
            else
            {
                m_matrix_rhs.diag(N - 1) = 2.0;
                m_matrix_lhs.diag(N - 1) = 2.0;
            }
        }

        // Infinite-dose donor: clamp every donor cell at its initial value.
        // Combined with the alpha override above, this makes the donor act as
        // a true Dirichlet reservoir at the donor / first-skin interface
        // regardless of D_donor or donor mesh density.
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
                if (i < N - 1)
                {
                    m_matrix_rhs.upper(i) = 0.0;
                    m_matrix_lhs.upper(i) = 0.0;
                }
            }
        }

        return true;
    }
}
