#ifndef SC_GEOMETRY_H
#define SC_GEOMETRY_H

#include "compartment.h"
#include "sink.h"

#include <vector>

namespace sc
{
    // Per-compartment uniform mesh, with cell size scaled by sqrt(D_i / D_min).
    // The smallest-D compartment uses cells of size 1/ss_per_um; higher-D
    // compartments get proportionally coarser cells, since their gradient
    // scale in the activity variable u = c/K is correspondingly wider. The
    // activity-FVM scheme is exact for piecewise-linear u, so refining at
    // compartment interfaces buys nothing in the bulk physics.
    class Geometry
    {
      public:
        Geometry() = default;

        // Builds the space-step vector and assigns geometry indices to each
        // compartment (and the sink, if non-null). Returns true on success.
        bool create(std::vector<Compartment>& compartments, int ss_per_um,
                    Sink* sink = nullptr);

        // Drops the half-open range [from_idx, to_idx) from the space-step vector.
        // Used after the donor compartment is removed mid-simulation.
        void remove(int from_idx, int to_idx);

        [[nodiscard]] const std::vector<double>& spaceSteps() const noexcept
        {
            return m_space_steps;
        }
        [[nodiscard]] int    size() const noexcept
        {
            return static_cast<int>(m_space_steps.size());
        }
        [[nodiscard]] double minSpaceStep() const noexcept { return m_min_space_step; }
        [[nodiscard]] double maxSpaceStep() const noexcept { return m_max_space_step; }
        [[nodiscard]] bool   valid() const noexcept { return m_valid; }

      private:
        std::vector<double> m_space_steps;
        double m_min_space_step = 1.0;
        double m_max_space_step = 1.0;
        bool   m_valid          = false;
    };
}

#endif  // SC_GEOMETRY_H
