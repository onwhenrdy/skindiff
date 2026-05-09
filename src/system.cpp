#include "system.h"

#include "algorithms.h"

#include <cassert>
#include <utility>

namespace sc
{
    namespace
    {
        // Convert cm^2 -> um^2.
        constexpr double cm2_to_um2(double cm2) noexcept { return cm2 * 1.0e8; }
        // Convert mg/ml -> mg/um^3 (1 ml = 1e12 um^3).
        constexpr double mg_per_ml_to_mg_per_um3(double v) noexcept { return v * 1.0e-12; }
        // Convert ml -> um^3.
        constexpr double ml_to_um3(double v) noexcept { return v * 1.0e12; }

        // Concentration `c` at cell idx, regardless of internal representation.
        // For c-formulation methods K_per_cell is all 1's, so this is the
        // identity. For activity methods (u-formulation) we recover c via
        // c = u * K.
        inline double cellConc(int idx, const std::vector<double>& state,
                               const std::vector<double>& K_per_cell) noexcept
        {
            return state[static_cast<std::size_t>(idx)] *
                   K_per_cell[static_cast<std::size_t>(idx)];
        }

        double integrateMass(const Compartment& comp, const Geometry& geometry,
                             const std::vector<double>& state,
                             const std::vector<double>& K_per_cell, double scale)
        {
            const auto& ss = geometry.spaceSteps();
            double mass = 0.0;
            for (int i = comp.geo_from; i <= comp.geo_to; ++i)
            {
                mass += cellConc(i, state, K_per_cell) * ss[static_cast<std::size_t>(i)];
            }
            return mass * scale * comp.area_um2;
        }

        double sinkMassValue(const Sink& sink, const Geometry& geometry,
                             const std::vector<double>& state,
                             const std::vector<double>& K_per_cell, double scale)
        {
            // The sink "cell" is virtual: the stored value (whether c or u) is
            // scaled by Vd/cell_vol (see initConcentrations) so that
            //     stored * K_sink * cell_vol = true mass in the PK compartment.
            // K_sink is 1 by convention, so cellConc(sink, ...) is the
            // Vd-scaled value and the integrand is the true mass.
            const auto idx  = sink.geo_from;
            const auto ss   = geometry.spaceSteps()[static_cast<std::size_t>(idx)];
            return cellConc(idx, state, K_per_cell) * ss * sink.area_um2 * scale;
        }

        // Sample concentration profile for a compartment in scaling units / ml.
        std::vector<double> sampleProfile(const Compartment& comp,
                                          const std::vector<double>& state,
                                          const std::vector<double>& K_per_cell, double scale)
        {
            const auto n = static_cast<std::size_t>(comp.geo_to - comp.geo_from + 1);
            std::vector<double> result(n);
            for (std::size_t k = 0; k < n; ++k)
            {
                const auto idx = comp.geo_from + static_cast<int>(k);
                result[k] = cellConc(idx, state, K_per_cell) * scale * 1.0e12;
            }
            return result;
        }

        // Cumulative mid-point depths (um) for the cells of a compartment, measured
        // from the top of the compartment.
        std::vector<double> compartmentDepths(const Compartment& comp, const Geometry& geometry)
        {
            const auto& ss = geometry.spaceSteps();
            std::vector<double> depths;
            depths.reserve(static_cast<std::size_t>(comp.geo_to - comp.geo_from + 1));

            double pos = 0.0;
            for (int i = comp.geo_from; i <= comp.geo_to; ++i)
            {
                const auto step = ss[static_cast<std::size_t>(i)];
                depths.push_back(pos + step / 2.0);
                pos += step;
            }
            return depths;
        }
    }

    System::System(Parameters parameters)
        : m_parameters(std::move(parameters))
        , m_sim_time(m_parameters.sys.simulation_time)
        , m_replace_after(m_parameters.vehicle.replace_after)
        , m_remove_at(m_parameters.vehicle.remove_at)
        , m_scale(scaleFactor(m_parameters.log.scaling))
    {
        const auto& v   = m_parameters.vehicle;
        const auto& sys = m_parameters.sys;
        const auto& sk  = m_parameters.sink;
        const auto& pk  = m_parameters.pk;

        m_matrix_builder.setMaxModule(sys.max_module);

        // Vehicle / donor compartment.
        const auto app_area_um2 = cm2_to_um2(v.app_area);
        Compartment donor{v.height, v.D, 1.0, app_area_um2, v.name};
        donor.c_init      = mg_per_ml_to_mg_per_um3(v.c_init);
        donor.finite_dose = v.finite_dose;
        m_compartments.push_back(std::move(donor));

        // Skin layers.
        for (const auto& l : m_parameters.layers)
        {
            Compartment c{l.height, l.D, l.K, app_area_um2 * l.cross_section, l.name};
            c.c_init = mg_per_ml_to_mg_per_um3(l.c_init);
            m_compartments.push_back(std::move(c));
        }

        // Sink.
        m_sink.name     = sk.name;
        m_sink.type     = pk.enabled ? Sink::Type::PK_Compartment : Sink::Type::Perfect_Sink;
        m_sink.area_um2 = app_area_um2 *
            (m_parameters.layers.empty() ? 1.0 : m_parameters.layers.back().cross_section);
        m_sink.Vd     = sk.Vd;
        m_sink.t_half = pk.thalf * 60.0;     // hours -> minutes
        m_sink.c_init = mg_per_ml_to_mg_per_um3(sk.c_init);

        m_geometry.setEta(sys.eta);
        buildGeometryAndMatrices();
        initConcentrations();
        initLoggers();
    }

    void System::buildGeometryAndMatrices()
    {
        m_geometry.create(m_parameters.sys.disc_method, m_compartments,
                          m_parameters.sys.resolution, &m_sink);
        m_matrix_builder.buildMatrix(m_compartments, m_geometry, &m_sink);
    }

    void System::initConcentrations()
    {
        const auto N = static_cast<std::size_t>(m_geometry.size());

        // Per-cell K, used to convert between the stored activity u and the
        // physical concentration c = u * K. The sink uses K = 1 (its activity
        // and concentration coincide).
        m_K_per_cell.assign(N, 1.0);
        for (const auto& comp : m_compartments)
        {
            for (int i = comp.geo_from; i <= comp.geo_to; ++i)
            {
                m_K_per_cell[static_cast<std::size_t>(i)] = comp.K;
            }
        }
        m_K_per_cell[static_cast<std::size_t>(m_sink.geo_from)] = 1.0;

        // Initial activity in each compartment cell: u = c_init / K.
        m_concentrations.assign(N, 0.0);
        for (const auto& comp : m_compartments)
        {
            const auto u_init = comp.c_init / comp.K;
            for (int i = comp.geo_from; i <= comp.geo_to; ++i)
            {
                m_concentrations[static_cast<std::size_t>(i)] = u_init;
            }
        }

        // Sink Vd-scaling: store cell value such that
        //     stored * cell_volume = c_init_sink * Vd  (the true PK-compartment mass).
        // K_sink = 1 so this is the same expression in both u- and c-spaces.
        const auto idx  = m_sink.geo_from;
        const auto ss   = m_geometry.spaceSteps()[static_cast<std::size_t>(idx)];
        const auto Vd_u = ml_to_um3(m_sink.Vd);
        m_concentrations[static_cast<std::size_t>(idx)] =
            m_sink.c_init * Vd_u / (ss * m_sink.area_um2);
    }

    void System::initLoggers()
    {
        const auto& v   = m_parameters.vehicle;
        const auto& log = m_parameters.log;

        m_mass_series.clear();
        m_cdp_series.clear();
        m_mass_series.resize(m_compartments.size());
        m_cdp_series.resize(m_compartments.size());

        for (std::size_t i = 0; i < m_compartments.size(); ++i)
        {
            const auto is_vehicle  = (i == 0);
            const auto mass_enabled =
                is_vehicle ? v.log_mass : m_parameters.layers[i - 1].log_mass;
            const auto cdp_enabled =
                is_vehicle ? v.log_cdp : m_parameters.layers[i - 1].log_cdp;

            m_mass_series[i].enabled      = mass_enabled;
            m_mass_series[i].log_interval = log.mass_log_interval;
            m_mass_series[i].reserve_for_total(m_sim_time);

            m_cdp_series[i].enabled      = cdp_enabled;
            m_cdp_series[i].log_interval = log.cdp_log_interval;
            m_cdp_series[i].depths_um    = compartmentDepths(m_compartments[i], m_geometry);
            m_cdp_series[i].reserve_for_total(m_sim_time);
        }

        m_sink_mass.enabled      = m_parameters.sink.log_mass;
        m_sink_mass.log_interval = log.mass_log_interval;
        m_sink_mass.reserve_for_total(m_sim_time);
    }

    void System::recordAt(double t)
    {
        for (std::size_t i = 0; i < m_compartments.size(); ++i)
        {
            if (m_mass_series[i].should_log(t))
            {
                m_mass_series[i].record(
                    t, integrateMass(m_compartments[i], m_geometry, m_concentrations,
                                     m_K_per_cell, m_scale));
            }
            if (m_cdp_series[i].should_log(t))
            {
                m_cdp_series[i].record(
                    t, sampleProfile(m_compartments[i], m_concentrations, m_K_per_cell, m_scale));
            }
        }
        if (m_sink_mass.should_log(t))
        {
            m_sink_mass.record(t, sinkMassValue(m_sink, m_geometry, m_concentrations,
                                                m_K_per_cell, m_scale));
        }
    }

    void System::replaceTopCompartment()
    {
        const auto& top  = m_compartments.front();
        const auto u_init = top.c_init / top.K;
        for (int i = top.geo_from; i <= top.geo_to; ++i)
        {
            m_concentrations[static_cast<std::size_t>(i)] = u_init;
        }
    }

    void System::removeTopCompartment()
    {
        const auto top      = m_compartments.front();
        const auto top_size = top.geo_to + 1;

        m_compartments.erase(m_compartments.begin());
        if (!m_mass_series.empty()) m_mass_series.erase(m_mass_series.begin());
        if (!m_cdp_series.empty())  m_cdp_series.erase(m_cdp_series.begin());

        m_geometry.remove(top.geo_from, top.geo_to + 1);

        m_concentrations.erase(m_concentrations.begin(),
                               m_concentrations.begin() + top_size);
        m_K_per_cell.erase(m_K_per_cell.begin(),
                           m_K_per_cell.begin() + top_size);

        for (auto& c : m_compartments)
        {
            c.geo_from -= top_size;
            c.geo_to   -= top_size;
        }
        m_sink.geo_from -= top_size;
        m_sink.geo_to   -= top_size;

        m_matrix_builder.buildMatrix(m_compartments, m_geometry, &m_sink);
    }

    System::Result System::run()
    {
        if (!initRun())
        {
            return Result::Failed;
        }

        auto n_ts       = m_matrix_builder.timesteps();
        auto rhs_matrix = m_matrix_builder.matrixRhs();
        auto lhs_matrix = m_matrix_builder.matrixLhs();

        bool vehicle_removed    = false;
        const bool must_replace = m_replace_after != 0;
        const bool must_remove  = m_remove_at != 0;

        recordAt(0.0);

        for (int t = 1; t <= m_sim_time; ++t)
        {
            if (testForStop(t))
            {
                return Result::Stopped;
            }
            progressCallback(t);

            for (int ts = 1; ts <= n_ts; ++ts)
            {
                rhs_matrix.inlineMultiply(m_concentrations);
                algorithm::thomasReUseIP(lhs_matrix, m_concentrations);
            }

            if (must_replace && !vehicle_removed && t > 1 && t % m_replace_after == 0)
            {
                replaceTopCompartment();
            }

            if (must_remove && t == m_remove_at)
            {
                vehicle_removed = true;
                removeTopCompartment();
                rhs_matrix = m_matrix_builder.matrixRhs();
                lhs_matrix = m_matrix_builder.matrixLhs();
                n_ts       = m_matrix_builder.timesteps();
            }

            recordAt(static_cast<double>(t));
        }

        if (!tearDownRun())
        {
            return Result::Failed;
        }
        return Result::Executed;
    }
}
