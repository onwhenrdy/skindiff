#include "compartmentlog2d.h"
#include "sink.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

#ifdef DSKIN_CORE_USE_ZIP
#include "zstr.h"
#endif

namespace sc
{
    CompartmentLog2D::CompartmentLog2D(MatrixBuilder::Method method, double app_area)
        : m_col_sep("\t")
        , m_filename("logger.dat")
        , m_column1_name("time")
        , m_column2_name("mass")
        , m_registered_sink(nullptr)
        , m_registered_compartment(nullptr)
        , m_auto_log_enabled(true)
        , m_enabled(true)
        , m_zip(false)
        , m_log_interval(1)
        , m_app_area(app_area)
        , m_mb_method(method)
    {
    }

    CompartmentLog2D::CompartmentLog2D(MatrixBuilder::Method method, double app_area,
                                       const std::string& name)
        : m_name(name)
        , m_filename("logger.dat")
        , m_column1_name("time")
        , m_column2_name("mass")
        , m_registered_sink(nullptr)
        , m_registered_compartment(nullptr)
        , m_auto_log_enabled(true)
        , m_enabled(true)
        , m_zip(false)
        , m_log_interval(1)
        , m_app_area(app_area)
        , m_mb_method(method)
    {
    }

    bool CompartmentLog2D::writeToFile() const
    {
#ifdef DSKIN_CORE_USE_ZIP
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
#endif
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

    const std::vector<double>& CompartmentLog2D::xs() const
    {
        return m_x;
    }

    const std::vector<double>& CompartmentLog2D::ys() const
    {
        return m_y;
    }

    int CompartmentLog2D::size() const
    {
        return m_x.size();
    }

    double CompartmentLog2D::x(int idx) const
    {
        return m_x[idx];
    }

    double CompartmentLog2D::y(int idx) const
    {
        return m_y[idx];
    }

    const std::string& CompartmentLog2D::name() const
    {
        return m_name;
    }

    void CompartmentLog2D::setName(const std::string& name)
    {
        m_name = name;
    }

    void sc::CompartmentLog2D::log(double x, double y)
    {
        m_x.push_back(x);
        m_y.push_back(y);
    }

    std::string CompartmentLog2D::columnSeparator() const
    {
        return m_col_sep;
    }

    void CompartmentLog2D::setColumnSeparator(const std::string& col_sep)
    {
        m_col_sep = col_sep;
    }

    std::string CompartmentLog2D::filename() const
    {
        return m_filename;
    }

    void CompartmentLog2D::setFilename(const std::string& filename)
    {
        m_filename = filename;
    }

    std::string CompartmentLog2D::column1Name() const
    {
        return m_column1_name;
    }

    std::string CompartmentLog2D::column2Name() const
    {
        return m_column2_name;
    }

    void CompartmentLog2D::setColumn2Name(const std::string& column2_name)
    {
        m_column2_name = column2_name;
    }

    void CompartmentLog2D::setColumnNames(const std::string& column1_name,
                                          const std::string& column2_name)
    {
        m_column1_name = column1_name;
        m_column2_name = column2_name;
    }

    Sink* CompartmentLog2D::registeredSink() const
    {
        return m_registered_sink;
    }

    void CompartmentLog2D::registerSink(Sink* registered_sink)
    {
        m_registered_sink        = registered_sink;
        m_registered_compartment = nullptr;
    }

    void CompartmentLog2D::log(double x_val, const Geometry& geometry,
                               const std::vector<double>& concentrations, double scale_fac)
    {
        if (!m_auto_log_enabled || (static_cast<int>(x_val) % m_log_interval != 0))
        {
            return;
        }

        if (m_registered_sink)
        {
            auto const A = (m_mb_method == MatrixBuilder::Method::DSkin_1_5)
                               ? m_app_area
                               : m_registered_sink->A();
            const auto idx  = m_registered_sink->geometryFromIdx();
            const auto conc = concentrations[idx];
            const auto ss   = geometry.spaceSteps()[idx];
            // heres: concentration
            const auto mass = conc * ss * A * scale_fac / m_registered_sink->Vd();
            this->log(x_val, mass);
        }
        else if (m_registered_compartment)
        {
            const auto idx_from = m_registered_compartment->geometryFromIdx();
            const auto idx_to   = m_registered_compartment->geometryToIdx();

            auto const A = (m_mb_method == MatrixBuilder::Method::DSkin_1_5)
                               ? m_app_area
                               : m_registered_compartment->A();
            auto mass = 0.0;
            for (int i = idx_from; i <= idx_to; ++i)
            {
                const auto ss   = geometry.spaceSteps()[i];
                const auto conc = concentrations[i];
                mass += conc * ss;
            }

            this->log(x_val, mass * scale_fac * A);
        }
        else
        {
            this->log(x_val, 0.0);
        }
    }

    bool CompartmentLog2D::autoLogEnabled() const
    {
        return m_auto_log_enabled;
    }

    void CompartmentLog2D::setAutoLogEnabled(bool value)
    {
        m_auto_log_enabled = value;
    }

    Compartment* CompartmentLog2D::registeredCompartment() const
    {
        return m_registered_compartment;
    }

    void CompartmentLog2D::registerCompartment(Compartment* registered_compartment)
    {
        m_registered_compartment = registered_compartment;
        m_registered_sink        = nullptr;
    }

    bool CompartmentLog2D::enabled() const
    {
        return m_enabled;
    }

    void CompartmentLog2D::setEnabled(bool enabled)
    {
        m_enabled = enabled;
    }

    bool CompartmentLog2D::zip() const
    {
        return m_zip;
    }

    void CompartmentLog2D::setZip(bool zip)
    {
        m_zip = zip;
    }

    int CompartmentLog2D::logInterval() const
    {
        return m_log_interval;
    }

    void CompartmentLog2D::setLogInterval(int log_interval)
    {
        assert(log_interval >= 1);

        const auto old_li = m_log_interval;
        m_log_interval    = log_interval;
        if (log_interval < old_li)
        {
            this->reserve();
        }
    }

    void CompartmentLog2D::setTimeHint(int time_hint)
    {
        assert(time_hint > 0);

        const auto old_th = m_time_hint;
        m_time_hint       = time_hint;
        if (time_hint > old_th)
        {
            this->reserve();
        }
    }

    void CompartmentLog2D::reserve()
    {
        // first entry is always logged at time = 0

        const auto new_cap = 1 + static_cast<int>(std::floor(m_time_hint / m_log_interval));
        this->m_x.reserve(new_cap);
        this->m_y.reserve(new_cap);
    }

    void CompartmentLog2D::setColumn1Name(const std::string& column1_name)
    {
        m_column1_name = column1_name;
    }

    std::ostream& operator<<(std::ostream& out, const CompartmentLog2D& logger)
    {
        out << logger.m_column1_name << logger.m_col_sep << logger.m_column2_name << "\n";

        const int size = logger.m_x.size();
        if (size > 0)
        {
            out.precision(std::numeric_limits<double>::digits10 + 1);
            for (int i = 0; i < size - 1; ++i)
            {
                out << logger.m_x[i] << logger.m_col_sep << logger.m_y[i] << "\n";
            }
            out << logger.m_x[size - 1] << logger.m_col_sep << logger.m_y[size - 1];
        }
        return out;
    }
}
