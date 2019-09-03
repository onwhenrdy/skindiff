#ifndef SC_COMPARTMENTLOG3D_H
#define SC_COMPARTMENTLOG3D_H

#include <map>
#include <ostream>
#include <string>
#include <vector>
#include "matrixbuilder.h"

namespace sc
{
    class Compartment;
}

namespace sc
{
    class CompartmentLog3D
    {
      public:
        enum class CPosition
        {
            Left,
            Center,
            Right
        };

        CompartmentLog3D();
        explicit CompartmentLog3D(const std::string& name);

        // in 1/um size
        void setStepSizes(const std::vector<double>& data);

        // time in minutes, data in mg/um^3
        void log(double time, const std::vector<double>& data);

        const std::vector<int>& times() const;
        std::vector<double> space() const;
        const std::vector<std::vector<double>>& data() const;

        bool writeToFile() const;
        friend std::ostream& operator<<(std::ostream& out, const CompartmentLog3D& logger);

        const std::string& name() const;
        void setName(const std::string& name);

        std::string columnSeparator() const;
        void setColumnSeparator(const std::string& col_sep);

        std::string filename() const;
        void setFilename(const std::string& filename);

        bool autoLogEnabled() const;
        void setAutoLogEnabled(bool value);

        Compartment* registeredCompartment() const;
        void registerCompartment(Compartment* registered_compartment);

        void log(double time, const std::vector<double>& concentrations, double scale_fac);

        CPosition concentrationPosition() const;
        void setConcentrationPosition(const CPosition& concentration_position);
        void setConcentrationPosition(MatrixBuilder::Method from_method);

        bool enabled() const;
        void setEnabled(bool enabled);

        bool zip() const;
        void setZip(bool zip);

        // in minutes
        int logInterval() const;
        void setLogInterval(int log_interval);

        // in minutes
        void setTimeHint(int time_hint);

    private:
        void reserve();

    private:
        std::string m_name;
        std::string m_col_sep;
        std::string m_filename;
        CPosition m_concentration_position;

        // payload
        std::vector<double> m_step_sizes;
        std::vector<int> m_times;
        std::vector<std::vector<double>> m_data;

        // auto logging stuff
        Compartment* m_registered_compartment;
        bool m_auto_log_enabled;

        // state
        bool m_enabled;
        bool m_zip;
        int m_log_interval;
        // size hint
        int m_time_hint;
    };

    std::ostream& operator<<(std::ostream& out, const CompartmentLog3D& logger);
}

#endif  // COMPARTMENTLOG3D_H
