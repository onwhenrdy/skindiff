#ifndef SC_COMPARTMENTLOG2D_H
#define SC_COMPARTMENTLOG2D_H

#include "geometry.h"
#include "matrixbuilder.h"

#include <string>
#include <vector>

namespace sc
{
    class Sink;
    class Compartment;
}

namespace sc
{
    class CompartmentLog2D
    {
      public:
        CompartmentLog2D(MatrixBuilder::Method method, double app_area);
        CompartmentLog2D(MatrixBuilder::Method method, double app_area, const std::string& name);
        void log(double x, double y);

        bool writeToFile() const;
        friend std::ostream& operator<<(std::ostream& out, const CompartmentLog2D& logger);

        const std::vector<double>& xs() const;
        const std::vector<double>& ys() const;
        int size() const;
        double x(int idx) const;
        double y(int idx) const;

        const std::string& name() const;
        void setName(const std::string& name);

        std::string columnSeparator() const;
        void setColumnSeparator(const std::string& col_sep);

        std::string filename() const;
        void setFilename(const std::string& filename);

        std::string column1Name() const;
        void setColumn1Name(const std::string& column1_name);

        std::string column2Name() const;
        void setColumn2Name(const std::string& column2_name);
        void setColumnNames(const std::string& column1_name, const std::string& column2_name);

        Sink* registeredSink() const;
        void registerSink(Sink* registered_sink);

        bool autoLogEnabled() const;
        void setAutoLogEnabled(bool value);

        Compartment* registeredCompartment() const;
        void registerCompartment(Compartment* registered_compartment);

        bool enabled() const;
        void setEnabled(bool enabled);

        bool zip() const;
        void setZip(bool zip);

        // in minutes
        int logInterval() const;
        void setLogInterval(int log_interval);

        // in minutes
        void setTimeHint(int time_hint);

        void log(double x_val, const Geometry& geometry, const std::vector<double>& concentrations,
                 double scale_fac);

      private:
        void reserve();

      private:
        std::vector<double> m_x;
        std::vector<double> m_y;
        std::string m_name;
        std::string m_col_sep;
        std::string m_filename;
        std::string m_column1_name;
        std::string m_column2_name;

        // auto logging stuff
        Sink* m_registered_sink;
        Compartment* m_registered_compartment;
        bool m_auto_log_enabled;

        // state
        bool m_enabled;
        bool m_zip;
        int m_log_interval;
        // size hint
        int m_time_hint;

        // state
        double m_app_area;
        MatrixBuilder::Method m_mb_method;
    };

    std::ostream& operator<<(std::ostream& out, const CompartmentLog2D& logger);
}
#endif  // COMPARTMENTLOG2D_H
