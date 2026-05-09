#ifndef SC_GEOMETRY_H
#define SC_GEOMETRY_H

#include "compartment.h"
#include "sink.h"

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace sc
{
    class Geometry
    {
      public:
        enum class DiscMethod
        {
            EQUI_DIST,  // Equidistant mesh: same dx everywhere.
            B_AND_K,    // Grid refinement by Babucke & Kloker, 2009: small
                        // cells at compartment interfaces, uniform 1 um in
                        // the bulk.
            GRADED      // Per-compartment uniform mesh with dx_i scaled by
                        // sqrt(D_i / D_min). High-D layers get coarser
                        // cells; the smallest-D layer gets cells of size
                        // 1/resolution.
        };

        Geometry() = default;

        // Builds the space-step vector and assigns geometry indices to each
        // compartment (and the sink, if non-null). Returns true on success.
        bool create(DiscMethod method, std::vector<Compartment>& compartments,
                    int ss_per_um, Sink* sink = nullptr);

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
        [[nodiscard]] DiscMethod discMethod() const noexcept { return m_disc_method; }
        [[nodiscard]] bool   valid() const noexcept { return m_valid; }

        [[nodiscard]] double eta() const noexcept { return m_eta; }
        void setEta(double eta) noexcept { m_eta = eta; }
        [[nodiscard]] double calculatedEta() const noexcept { return m_calculated_eta; }

      private:
        void findOptTransition(int& n, double& x, int& a, double& delta_x, double err);
        double findOptimalX(double start_x, int n, double a, double err) const;
        static double powerSeriesDoubleLastElement(int n, double x) noexcept;

        std::vector<double> m_space_steps;
        double     m_min_space_step = 1.0;
        double     m_max_space_step = 1.0;
        DiscMethod m_disc_method    = DiscMethod::EQUI_DIST;
        bool       m_valid          = false;
        double     m_eta            = 0.6;
        double     m_calculated_eta = 0.0;
    };

    [[nodiscard]] std::string_view toString(Geometry::DiscMethod method) noexcept;
    [[nodiscard]] std::optional<Geometry::DiscMethod>
        discMethodFromString(std::string_view str) noexcept;
}

#endif  // SC_GEOMETRY_H
