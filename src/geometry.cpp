#include "geometry.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <limits>

namespace sc
{
    bool Geometry::create(Geometry::DiscMethod method, std::vector<Compartment>& compartments,
                          int ss_per_um, Sink* sink)
    {
        m_disc_method = method;
        const int c_size = static_cast<int>(compartments.size());
        assert(c_size > 0);

        m_space_steps.clear();
        m_min_space_step = 1.0;
        m_max_space_step = 1.0;

        if (method == DiscMethod::EQUI_DIST || ss_per_um == 1 || c_size == 1)
        {
            const auto ss = 1.0 / ss_per_um;
            int counter = 0;
            int total_um = 0;
            for (auto& c : compartments)
            {
                assert(c.height_um > 0);
                const auto start_idx = counter;
                counter += c.height_um * ss_per_um - 1;
                c.geo_from = start_idx;
                c.geo_to   = counter;
                ++counter;
                total_um += c.height_um;
            }

            int size = total_um * ss_per_um;
            if (sink)
            {
                ++size;
                sink->geo_from = counter;
                sink->geo_to   = counter;
            }

            m_valid = (size > 0);
            if (!m_valid) return false;

            m_min_space_step = ss;
            m_max_space_step = ss;
            m_space_steps    = std::vector<double>(size, ss);
            return true;
        }

        if (method == DiscMethod::B_AND_K)
        {
            constexpr double eps = 1.0e-13;
            m_calculated_eta     = m_eta;
            int  n_trans_eles    = 1;
            int  n_trans_size    = 1;
            double ss_boundary   = 1.0 / ss_per_um;
            findOptTransition(n_trans_eles, m_calculated_eta, n_trans_size, ss_boundary, eps);
            assert(n_trans_eles > 0);

            std::vector<double> trans_vec;
            trans_vec.reserve(static_cast<std::size_t>(n_trans_eles) * 2);
            double ss = 1.0;
            for (int i = 0; i < n_trans_eles - 1; ++i)
            {
                ss *= m_calculated_eta;
                trans_vec.push_back(ss);
            }
            trans_vec.push_back(ss);
            for (int i = n_trans_eles - 1; i >= 0; --i)
            {
                trans_vec.push_back(trans_vec[static_cast<std::size_t>(i)]);
            }

            int counter = 0;
            int c_carry = 0;
            const auto t_vec_size = static_cast<int>(trans_vec.size());
            for (int i = 0; i < c_size; ++i)
            {
                auto& c = compartments[static_cast<std::size_t>(i)];
                const auto start_idx = counter;
                const auto size      = c.height_um -
                    ((i == 0 || i == c_size - 1) ? n_trans_size : n_trans_size * 2);
                assert(size >= 0);
                for (int j = 0; j < size; ++j)
                {
                    m_space_steps.push_back(1.0);
                    ++counter;
                }
                if (i < c_size - 1)
                {
                    counter += t_vec_size / 2;
                    for (auto val : trans_vec)
                    {
                        m_space_steps.push_back(val);
                    }
                }
                counter += c_carry;
                c.geo_from = start_idx;
                c.geo_to   = counter - 1;
                c_carry    = t_vec_size / 2;
            }

            if (sink)
            {
                m_space_steps.push_back(1.0);
                sink->geo_from = counter;
                sink->geo_to   = counter;
            }

            m_valid          = true;
            m_max_space_step = 1.0;
            m_min_space_step = ss_boundary;
            return true;
        }

        if (method == DiscMethod::GRADED)
        {
            // Per-compartment uniform mesh, with dx scaled by sqrt(D_i / D_min).
            // Smallest-D compartment uses cells of size 1/ss_per_um; higher-D
            // compartments get proportionally coarser cells, since a larger D
            // gives a wider gradient scale in the activity variable u and we
            // don't need to resolve the bulk as finely. Cell count per
            // compartment is rounded so the cells fit the compartment exactly.
            assert(ss_per_um > 0);
            double D_min = compartments.front().D;
            for (const auto& c : compartments)
            {
                if (c.D > 0.0 && c.D < D_min) D_min = c.D;
            }
            assert(D_min > 0.0);

            const double dx_min = 1.0 / ss_per_um;

            int counter = 0;
            for (auto& c : compartments)
            {
                assert(c.height_um > 0);
                const auto start_idx = counter;

                const double dx_target = dx_min * std::sqrt(std::max(c.D, D_min) / D_min);
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

        m_valid = false;
        return false;
    }

    void Geometry::remove(int from_idx, int to_idx)
    {
        m_space_steps.erase(m_space_steps.begin() + from_idx, m_space_steps.begin() + to_idx);
        if (m_space_steps.empty()) return;

        const auto mm = std::minmax_element(m_space_steps.begin(), m_space_steps.end());
        m_min_space_step = *mm.first;
        m_max_space_step = *mm.second;
    }

    void Geometry::findOptTransition(int& n, double& x, int& a, double& delta_x, double err)
    {
        const double start_x = x;
        n = static_cast<int>(std::ceil(std::log10(delta_x) / std::log10(start_x)));
        a = static_cast<int>(std::ceil(powerSeriesDoubleLastElement(n, start_x)));
        x = findOptimalX(start_x, n, a, err);
        int i = 0;
        while (std::pow(x, n - 1) > delta_x && i < 10)
        {
            ++i;
            ++n;
            a = static_cast<int>(std::ceil(powerSeriesDoubleLastElement(n, start_x)));
            x = findOptimalX(start_x, n, a, err);
        }
        delta_x = std::pow(x, n - 1);
    }

    double Geometry::findOptimalX(double start_x, int n, double a, double err) const
    {
        auto x     = start_x;
        auto old_x = x + 2 * err;

        const auto dx = std::numeric_limits<double>::epsilon();
        while (std::abs(old_x - x) > err)
        {
            old_x          = x;
            const auto f_x  = powerSeriesDoubleLastElement(n, x) - a;
            const auto f_dx = powerSeriesDoubleLastElement(n, x + dx) - a;
            x -= f_x * dx / (f_dx - f_x);
        }
        return x;
    }

    double Geometry::powerSeriesDoubleLastElement(int n, double x) noexcept
    {
        double sum = 0.0;
        double last_comp = 1.0;
        for (int i = 1; i < n; ++i)
        {
            last_comp *= x;
            sum += last_comp;
        }
        sum += last_comp;
        return sum;
    }

    std::string_view toString(Geometry::DiscMethod method) noexcept
    {
        switch (method)
        {
            case Geometry::DiscMethod::EQUI_DIST: return "EQUIDIST";
            case Geometry::DiscMethod::B_AND_K:   return "BK";
            case Geometry::DiscMethod::GRADED:    return "GRADED";
        }
        return "EQUIDIST";
    }

    std::optional<Geometry::DiscMethod> discMethodFromString(std::string_view str) noexcept
    {
        std::string upper(str);
        std::transform(upper.begin(), upper.end(), upper.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        if (upper == "EQUIDIST") return Geometry::DiscMethod::EQUI_DIST;
        if (upper == "BK")       return Geometry::DiscMethod::B_AND_K;
        if (upper == "GRADED")   return Geometry::DiscMethod::GRADED;
        return std::nullopt;
    }
}
