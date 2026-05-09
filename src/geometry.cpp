#include "geometry.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

namespace sc
{
    bool Geometry::create(std::vector<Compartment>& compartments, int ss_per_um, Sink* sink)
    {
        const int c_size = static_cast<int>(compartments.size());
        assert(c_size > 0);
        assert(ss_per_um > 0);

        m_space_steps.clear();

        double D_min = std::numeric_limits<double>::infinity();
        for (const auto& c : compartments)
        {
            if (c.D > 0.0 && c.D < D_min) D_min = c.D;
        }
        assert(std::isfinite(D_min) && D_min > 0.0);

        const double dx_min = 1.0 / ss_per_um;

        int counter = 0;
        for (auto& c : compartments)
        {
            assert(c.height_um > 0);
            const auto start_idx = counter;

            const double dx_target = dx_min * std::sqrt(c.D / D_min);
            int n_cells = static_cast<int>(std::round(c.height_um / dx_target));
            if (n_cells < 1) n_cells = 1;
            const double actual_dx =
                static_cast<double>(c.height_um) / static_cast<double>(n_cells);

            for (int j = 0; j < n_cells; ++j)
            {
                m_space_steps.push_back(actual_dx);
                ++counter;
            }
            c.geo_from = start_idx;
            c.geo_to   = counter - 1;
        }

        if (sink)
        {
            m_space_steps.push_back(dx_min);
            sink->geo_from = counter;
            sink->geo_to   = counter;
        }

        const auto mm = std::minmax_element(m_space_steps.begin(), m_space_steps.end());
        m_min_space_step = *mm.first;
        m_max_space_step = *mm.second;
        m_valid          = true;
        return true;
    }

    void Geometry::remove(int from_idx, int to_idx)
    {
        m_space_steps.erase(m_space_steps.begin() + from_idx, m_space_steps.begin() + to_idx);
        if (m_space_steps.empty()) return;

        const auto mm = std::minmax_element(m_space_steps.begin(), m_space_steps.end());
        m_min_space_step = *mm.first;
        m_max_space_step = *mm.second;
    }
}
