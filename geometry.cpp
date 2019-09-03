#include "geometry.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <numeric>

#include "helper.h"

namespace sc
{
    Geometry::Geometry()
        : m_min_space_step(1.0)
        , m_max_space_step(1.0)
        , m_disc_method(Geometry::DiscMethod::EQUI_DIST)
        , m_valid(false)
        , m_eta(0.6)
        , m_calcuated_eta(0.0)
    {
    }

    bool Geometry::create(Geometry::DiscMethod method, std::vector<Compartment>& compartments,
                          int ss_per_um, Sink* sink)
    {
        m_disc_method = method;
        const int c_size = compartments.size();
        assert(c_size > 0);

        m_space_steps.clear();
        m_min_space_step = 1.0;
        m_max_space_step = 1.0;

        if (method == Geometry::DiscMethod::EQUI_DIST || ss_per_um == 1 || compartments.size() == 1)
        {
            auto ss     = 1.0 / ss_per_um;
            int size    = 0;
            int counter = 0;
            for (auto& c : compartments)
            {
                const auto comp_size = c.size();
                assert(comp_size > 0);

                const auto start_idx = counter;
                counter += comp_size * ss_per_um - 1;
                c.setGeometryIdx(start_idx, counter);
                counter++;
                size += comp_size;
            }

            size *= ss_per_um;
            if (sink)
            {
                size++;
                sink->setGeometryIdx(counter, counter);
            }

            m_valid = (size > 0);
            if (m_valid)
            {
                m_min_space_step = ss;
                m_max_space_step = ss;
                m_space_steps    = std::vector<double>(size, ss);
            }
            else
            {
                return false;
            }

            return true;
        }

        if (method == Geometry::DiscMethod::B_AND_K)
        {
            // find transition parameters
            static const auto eps = 1.0E-13;
            m_calcuated_eta       = m_eta;  // transition scaling factor
            auto n_trans_eles     = 1;      // # transition elements
            auto n_trans_size     = 1;      // transition size in um
            auto ss_boundary      = 1.0 / ss_per_um;
            this->findOptTransition(n_trans_eles, m_calcuated_eta, n_trans_size, ss_boundary, eps);
            assert(n_trans_eles > 0);

            // build transition vector (symetric)
            std::vector<double> trans_vec;
            trans_vec.reserve(n_trans_eles * 2);
            auto ss = 1.0;
            for (auto i = 0; i < n_trans_eles - 1; ++i)
            {
                ss *= m_calcuated_eta;
                trans_vec.push_back(ss);
            }
            trans_vec.push_back(ss);
            for (auto i = n_trans_eles - 1; i >= 0; --i)
            {
                trans_vec.push_back(trans_vec[i]);
            }

            // build the ss-vector
            int counter           = 0;
            int c_carry           = 0;
            const auto t_vec_size = trans_vec.size();
            for (auto i = 0; i < c_size; ++i)
            {
                auto& c              = compartments[i];
                const auto start_idx = counter;
                const auto size =
                    c.size() - ((i == 0 || i == c_size - 1) ? n_trans_size : n_trans_size * 2);
                assert(size >= 0);
                for (auto j = 0; j < size; ++j)
                {
                    m_space_steps.push_back(1.0);
                    counter++;
                }
                // add full transition
                if (i < c_size - 1)
                {
                    counter += t_vec_size / 2;
                    for (auto val : trans_vec)
                    {
                        m_space_steps.push_back(val);
                    }
                }
                counter += c_carry;
                const auto end_idx = counter - 1;
                c_carry            = t_vec_size / 2;
                c.setGeometryIdx(start_idx, end_idx);
            }

            if (sink)
            {
                m_space_steps.push_back(1.0);
                sink->setGeometryIdx(counter, counter);
            }

            m_valid          = true;
            m_max_space_step = 1.0;
            m_min_space_step = ss_boundary;
            return true;
        }

        m_valid = false;
        return false;
    }

    const std::vector<double>& Geometry::spaceSteps() const
    {
        return m_space_steps;
    }

    int Geometry::size() const
    {
        return m_space_steps.size();
    }

    std::string Geometry::spaceStepsR(const std::string& var_name) const
    {
        return toRVector(m_space_steps, var_name);
    }

    double Geometry::minSpaceStep() const
    {
        return m_min_space_step;
    }

    double Geometry::maxSpaceStep() const
    {
        return m_max_space_step;
    }

    Geometry::DiscMethod Geometry::discMethod() const
    {
        return m_disc_method;
    }

    bool Geometry::valid() const
    {
        return m_valid;
    }

    void Geometry::remove(int from_idx, int to_idx)
    {
        m_space_steps.erase(m_space_steps.begin() + from_idx, m_space_steps.begin() + to_idx);

        const auto min_max_it = std::minmax_element(m_space_steps.begin(), m_space_steps.end());
        m_max_space_step      = *min_max_it.second;
        m_min_space_step      = *min_max_it.first;
    }

    void Geometry::findOptTransition(int& n, double& x, int& a, double& delta_x, double err)
    {
        // two equal elements at the boundary
        double start_x = x;
        n              = static_cast<int>(std::ceil(log10(delta_x) / log10(start_x)));
        a              = static_cast<int>(std::ceil(powerSeriesDoubleLastElement(n, start_x)));
        x              = findOptimalX(start_x, n, a, err);
        int i          = 0;
        while (pow(x, n - 1) > delta_x && i < 10)
        {
            i += 1;
            n += 1;
            a = static_cast<int>(std::ceil(powerSeriesDoubleLastElement(n, start_x)));
            x = this->findOptimalX(start_x, n, a, err);
        }
        delta_x = std::pow(x, n - 1);
    }

    double Geometry::findOptimalX(double start_x, int n, double a, double err)
    {
        auto x     = start_x;
        auto old_x = x + 2 * err;

        const auto dx = std::numeric_limits<double>::epsilon();
        // find with newton algorithm
        while (std::abs(old_x - x) > err)
        {
            old_x           = x;
            const auto f_x  = this->powerSeriesDoubleLastElement(n, x) - a;
            const auto f_dx = this->powerSeriesDoubleLastElement(n, x + dx) - a;
            x -= f_x * dx / (f_dx - f_x);
        }
        return x;
    }

    double Geometry::powerSeriesDoubleLastElement(int n, double x)
    {
        auto sum       = 0.0;
        auto last_comp = 1.0;
        for (int i = 1; i < n; i++)
        {
            last_comp *= x;
            sum += last_comp;
        }
        sum += last_comp;
        return sum;
    }

    double Geometry::calculatedEta() const
    {
        return m_calcuated_eta;
    }

    double Geometry::eta() const
    {
        return m_eta;
    }

    void Geometry::setEta(double eta)
    {
        m_eta = eta;
    }

    std::string toString(Geometry::DiscMethod method)
    {
        switch (method)
        {
            case Geometry::DiscMethod::EQUI_DIST:
                return "EQUIDIST";
                break;
            case Geometry::DiscMethod::B_AND_K:
                return "BK";
                break;
            default:
                return "unknown";
                break;
        }
    }

    Geometry::DiscMethod discMethodFromString(const std::string& str, bool* ok)
    {
        auto toUpr = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
                return static_cast<unsigned char>(std::toupper(c));
            });
            return s;
        };

        if (ok)
        {
            *ok = true;
        }

        auto upp_str = toUpr(str);
        if (upp_str == "EQUIDIST")
        {
            return Geometry::DiscMethod::EQUI_DIST;
        }

        if (upp_str == "BK")
        {
            return Geometry::DiscMethod::B_AND_K;
        }

        if (ok)
        {
            *ok = false;
        }

        return Geometry::DiscMethod::EQUI_DIST;
    }
}
