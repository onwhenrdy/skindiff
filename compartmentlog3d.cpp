#include "compartmentlog3d.h"

#include "compartment.h"
#include "zstr.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <limits>

namespace sc
{
    CompartmentLog3D::CompartmentLog3D()
        : m_col_sep("\t")
        , m_filename("unknown.dat")
        , m_concentration_position(CompartmentLog3D::CPosition::Left)
        , m_registered_compartment(nullptr)
        , m_auto_log_enabled(true)
        , m_enabled(false)
        , m_zip(true)
        , m_log_interval(1)
    {
    }

    CompartmentLog3D::CompartmentLog3D(const std::string& name)
        : m_name(name)
        , m_col_sep("\t")
        , m_filename("unknown.dat")
        , m_concentration_position(CompartmentLog3D::CPosition::Left)
        , m_registered_compartment(nullptr)
        , m_auto_log_enabled(true)
        , m_enabled(false)
        , m_zip(true)
        , m_log_interval(1)
    {
    }

    void CompartmentLog3D::setStepSizes(const std::vector<double>& data)
    {
        m_step_sizes = data;
    }

    void CompartmentLog3D::log(double time, const std::vector<double>& data)
    {
        m_times.push_back(time);
        m_data.push_back(data);
    }

    const std::vector<int>& CompartmentLog3D::times() const
    {
        return m_times;
    }

    std::vector<double> CompartmentLog3D::space() const
    {
        std::vector<double> result;
        result.reserve(m_step_sizes.size());

        auto x_pos   = 0.0;
        auto pre_inc = 0.0;
        auto pos_inc = 0.0;
        for (auto val : m_step_sizes)
        {
            if (m_concentration_position == CompartmentLog3D::CPosition::Center)
            {
                pre_inc = val / 2.0;
                pos_inc = pre_inc;
            }
            else
            {
                pre_inc                                                                   = 0.0;
                pos_inc                                                                   = 0.0;
                (m_concentration_position == CompartmentLog3D::CPosition::Left ? pos_inc
                                                                               : pre_inc) = val;
            }

            x_pos += pre_inc;
            result.push_back(x_pos);
            x_pos += pos_inc;
        }

        return result;
    }

    const std::vector<std::vector<double>>& CompartmentLog3D::data() const
    {
        return m_data;
        ;
    }

    bool CompartmentLog3D::writeToFile() const
    {
        if (m_zip)
        {
            zstr::ofstream file(m_filename + ".gz");
            if (!file)
            {
                return false;
            }

            file << *this;
        }
        else
        {
            std::ofstream file(m_filename);
            if (!file)
            {
                return false;
            }

            file << *this;
        }

        return true;
    }

    const std::string& CompartmentLog3D::name() const
    {
        return m_name;
    }

    void CompartmentLog3D::setName(const std::string& name)
    {
        m_name = name;
    }

    std::string CompartmentLog3D::columnSeparator() const
    {
        return m_col_sep;
    }

    void CompartmentLog3D::setColumnSeparator(const std::string& col_sep)
    {
        m_col_sep = col_sep;
    }

    std::string CompartmentLog3D::filename() const
    {
        return m_filename;
    }

    void CompartmentLog3D::setFilename(const std::string& filename)
    {
        m_filename = filename;
    }

    bool CompartmentLog3D::autoLogEnabled() const
    {
        return m_auto_log_enabled;
    }

    void CompartmentLog3D::setAutoLogEnabled(bool value)
    {
        m_auto_log_enabled = value;
    }

    Compartment* CompartmentLog3D::registeredCompartment() const
    {
        return m_registered_compartment;
    }

    void CompartmentLog3D::registerCompartment(Compartment* registered_compartment)
    {
        m_registered_compartment = registered_compartment;
    }

    void CompartmentLog3D::log(double time, const std::vector<double>& concentrations,
                               double scale_fac)
    {
        if (!m_auto_log_enabled || (static_cast<int>(time) % m_log_interval != 0))
        {
            return;
        }

        auto data = std::vector<double>(m_step_sizes.size(), 0.0);
        if (m_registered_compartment)
        {
            const auto idx_from = m_registered_compartment->geometryFromIdx();
            const auto idx_to   = m_registered_compartment->geometryToIdx();
            assert(data.size() == (idx_to - idx_from) + 1);

            int idx = 0;
            for (int i = idx_from; i <= idx_to; ++i)
            {
                data[idx] = concentrations[i] * scale_fac;
                idx++;
            }
        }
        this->log(time, data);
    }

    CompartmentLog3D::CPosition CompartmentLog3D::concentrationPosition() const
    {
        return m_concentration_position;
    }

    void CompartmentLog3D::setConcentrationPosition(const CPosition& concentration_position)
    {
        m_concentration_position = concentration_position;
    }

    void CompartmentLog3D::setConcentrationPosition(MatrixBuilder::Method from_method)
    {
        switch (from_method)
        {
            case MatrixBuilder::Method::DSkin_1_3:
                m_concentration_position = CompartmentLog3D::CPosition::Center;
                break;
            case MatrixBuilder::Method::DSkin_1_4:
                m_concentration_position = CompartmentLog3D::CPosition::Center;
                break;
            case MatrixBuilder::Method::DSkin_1_5:
                m_concentration_position = CompartmentLog3D::CPosition::Center;
                break;
            default:
                assert(false);
                break;
        }
    }

    bool CompartmentLog3D::enabled() const
    {
        return m_enabled;
    }

    void CompartmentLog3D::setEnabled(bool enabled)
    {
        m_enabled = enabled;
    }

    bool CompartmentLog3D::zip() const
    {
        return m_zip;
    }

    void CompartmentLog3D::setZip(bool zip)
    {
        m_zip = zip;
    }

    int CompartmentLog3D::logInterval() const
    {
        return m_log_interval;
    }

    void CompartmentLog3D::setLogInterval(int log_interval)
    {
        assert(log_interval >= 1);

        const auto old_li = m_log_interval;
        m_log_interval    = log_interval;
        if (log_interval < old_li)
        {
            this->reserve();
        }
    }

    void CompartmentLog3D::setTimeHint(int time_hint)
    {
        assert(time_hint > 0);

        const auto old_th = m_time_hint;
        m_time_hint       = time_hint;
        if (time_hint > old_th)
        {
            this->reserve();
        }
    }

    void CompartmentLog3D::reserve()
    {
        // first entry is always logged at time = 0

        const auto new_cap = 1 + static_cast<int>(std::floor(m_time_hint / m_log_interval));
        this->m_times.reserve(new_cap);
        this->m_data.reserve(new_cap);
    }

    std::ostream& operator<<(std::ostream& out, const CompartmentLog3D& logger)
    {
        // format is
        // 0    SEP     X_POS_1     SEP     X_POS_2 ... \n
        // TIME_1   SEP     C_1         SEP     C_2     ... \n
        // ...

        out.precision(std::numeric_limits<double>::digits10 + 1);
        const auto sep = logger.columnSeparator();
        // x-vals
        out << 0;
        auto x_pos     = 0.0;
        const auto c_p = logger.m_concentration_position;
        auto pre_inc   = 0.0;
        auto pos_inc   = 0.0;
        for (auto val : logger.m_step_sizes)
        {
            if (c_p == CompartmentLog3D::CPosition::Center)
            {
                pre_inc = val / 2.0;
                pos_inc = pre_inc;
            }
            else
            {
                pre_inc                                                        = 0.0;
                pos_inc                                                        = 0.0;
                (c_p == CompartmentLog3D::CPosition::Left ? pos_inc : pre_inc) = val;
            }

            x_pos += pre_inc;
            out << sep << x_pos;
            x_pos += pos_inc;
        }
        out << "\n";

        // conc-vals
        const auto data_size = logger.m_data.size();
        int i                = 1;
        for (const auto& data : logger.m_data)
        {
            out << logger.m_times[i - 1];
            for (auto val : data)
            {
                out << sep << val;
            }

            if (i != data_size)
            {
                out << "\n";
            }
            i++;
        }

        return out;
    }
}
