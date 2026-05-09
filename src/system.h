#ifndef SC_SYSTEM_H
#define SC_SYSTEM_H

#include "compartment.h"
#include "geometry.h"
#include "logger.h"
#include "matrixbuilder.h"
#include "parameter.h"
#include "sink.h"

#include <vector>

namespace sc
{
    // Owns the discretized stack and the time-series loggers and runs the
    // Crank-Nicolson time-stepping loop.
    class System
    {
      public:
        enum class Result
        {
            Executed,
            Failed,
            Stopped
        };

        explicit System(Parameters parameters);
        virtual ~System() = default;

        Result run();

        [[nodiscard]] const Parameters& parameters() const noexcept { return m_parameters; }
        [[nodiscard]] const Geometry&   geometry()   const noexcept { return m_geometry; }
        [[nodiscard]] const std::vector<Compartment>& compartments() const noexcept
        {
            return m_compartments;
        }
        [[nodiscard]] const Sink& sink() const noexcept { return m_sink; }
        [[nodiscard]] const std::vector<double>& concentrations() const noexcept
        {
            return m_concentrations;
        }

        [[nodiscard]] const std::vector<MassSeries>& compartmentMass() const noexcept
        {
            return m_mass_series;
        }
        [[nodiscard]] const MassSeries& sinkMass() const noexcept { return m_sink_mass; }
        [[nodiscard]] const std::vector<CdpSeries>& cdp() const noexcept { return m_cdp_series; }

      protected:
        // Hooks for derived classes (e.g. R bindings) to inject progress / cancellation.
        virtual bool initRun()                       { return true; }
        virtual bool tearDownRun()                   { return true; }
        virtual void progressCallback(int /*t*/)     {}
        virtual bool testForStop(int /*t*/)          { return false; }

      private:
        void buildGeometryAndMatrices();
        void initConcentrations();
        void initLoggers();
        void recordAt(double t);
        void replaceTopCompartment();
        void removeTopCompartment();

        Parameters               m_parameters;
        std::vector<Compartment> m_compartments;
        // Internal state: the activity u = c/K, in units of mg/um^3.
        // Convert to physical concentration c via c_i = m_concentrations[i] *
        // m_K_per_cell[i]. The sink cell follows the Vd-scaling convention
        // (cell mass = c_sink_pk * Vd).
        std::vector<double>      m_concentrations;
        // Per-cell partition coefficient K (= K_layer per skin cell, K = 1
        // in the sink). Multiply m_concentrations[i] by m_K_per_cell[i] to
        // recover the physical concentration.
        std::vector<double>      m_K_per_cell;
        Sink                     m_sink;
        Geometry                 m_geometry;
        MatrixBuilder            m_matrix_builder;

        std::vector<MassSeries> m_mass_series;  // one per active compartment (incl. donor)
        MassSeries              m_sink_mass;
        std::vector<CdpSeries>  m_cdp_series;   // one per active compartment

        int    m_sim_time      = 1;
        int    m_replace_after = 0;
        int    m_remove_at     = 0;
        double m_scale         = 1.0;
    };
}

#endif  // SC_SYSTEM_H
