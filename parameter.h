#ifndef SC_PARAMETER_H
#define SC_PARAMETER_H

#include "geometry.h"
#include "matrixbuilder.h"
#include <string>

namespace sc
{
    // System parameter
    class SystemParameter
    {
      public:

        SystemParameter();
        bool isValid(std::string& error_string) const;
        std::string overviewString() const;

        Geometry::DiscMethod discMethod() const;
        void setDiscMethod(const Geometry::DiscMethod& disc_method);

        // in 1/x um
        int resolution() const;
        void setResolution(int resolution);

        // in minutes
        int simulationTime() const;
        void setSimulationTime(int simulation_time);

        double maxModule() const;
        void setMaxModule(double max_module);

        MatrixBuilder::Method matrixBuilderMethod() const;
        void setMatrixBuilderMethod(const MatrixBuilder::Method& matrix_builder_method);

        double eta() const;
        void setEta(double eta);

    private:
        Geometry::DiscMethod m_disc_method;
        MatrixBuilder::Method m_matrix_builder_method;
        int m_resolution;
        double m_max_module;
        double m_eta;
        int m_simulation_time;
    };

    // Pk parameter
    class PKParameter
    {
      public:
        PKParameter();

        bool isValid(std::string& error_string) const;
        std::string overviewString() const;

        bool enabled() const;
        void setEnabled(bool enabled);

        // in h
        double thalf() const;
        void setThalf(double thalf);

      private:
        bool m_enabled;
        double m_thalf;
    };

    class SinkParameter
    {
      public:
        SinkParameter();

        bool isValid(std::string& error_string) const;
        std::string overviewString() const;

        bool log() const;
        void setLog(bool log);

        const std::string& name() const;
        void setName(const std::string& name);

        // in ml
        double Vd() const;
        void setVd(double Vd);

        // in mg/ml
        double cInit() const;
        void setCInit(double value);

      private:
        bool m_log;
        std::string m_name;
        double m_Vd;
        double m_c_init;
    };

    class VehicleParameter
    {
      public:
        VehicleParameter();

        bool isValid(std::string& error_string) const;
        std::string overviewString() const;

        bool log() const;
        void setLog(bool log);

        const std::string& name() const;
        void setName(const std::string& name);

        // in mg/ml
        double cInit() const;
        void setCInit(double value);

        // in cm^2
        double appArea() const;
        void setAppArea(double app_area);

        // in um
        int height() const;
        void setHeight(int height);

        // in min
        int replaceAfter() const;
        void setReplaceAfter(int replace_after);
        bool replace() const;

        // in min
        int removeAt() const;
        void setRemoveAt(int remove_at);
        bool remove() const;

        // in um^2/min
        double D() const;
        void setD(double D);

        bool finiteDose() const;
        void setFiniteDose(bool finite_dose);

        bool logCDP() const;
        void setLogCDP(bool log_cdp);

      private:
        bool m_log;
        bool m_log_cdp;
        std::string m_name;
        double m_c_init;
        double m_app_area;
        double m_D;
        int m_height;
        int m_replace_after;
        int m_remove_at;
        bool m_finite_dose;
    };

    class LayerParameter
    {
      public:
        LayerParameter();

        bool isValid(std::string& error_string) const;
        std::string overviewString() const;

        bool log() const;
        void setLog(bool log);

        const std::string& name() const;
        void setName(const std::string& name);

        // in mg/ml
        double cInit() const;
        void setCInit(double value);

        // in um
        int height() const;
        void setHeight(int height);

        // in um^2/min
        double D() const;
        void setD(double D);

        // in K_Layer/Vehicle
        double K() const;
        void setK(double K);

        // in ]0..1]
        double crossSection() const;
        void setCrossSection(double cross_section);

        bool logCDP() const;
        void setLogCDP(bool log_cdp);

      private:
        bool m_log;
        bool m_log_cdp;
        std::string m_name;
        double m_c_init;
        double m_D;
        double m_K;
        double m_cross_section;
        int m_height;
    };

    class LogParameter
    {
      public:
        enum class Scaling
        {
            MG,
            UG,
            NG
        };

        LogParameter();
        bool isValid(std::string& error_string) const;
        std::string overviewString() const;

        bool gzipCDP() const;
        void setGzipCDP(bool gzip_CDP);

        bool gzipMass() const;
        void setGzipMass(bool gzip_mass);

        // in minutes
        int massLogInterval() const;
        void setMassLogInterval(int mass_log_interval);

        // in minutes
        int CDPLogInterval() const;
        void setCDPLogInterval(int cdp_log_interval);

        const std::string& massFilePostfix() const;
        void setMassFilePostfix(const std::string& mass_file_postfix);

        const std::string& CDPFilePostfix() const;
        void setCDPFilePostfix(const std::string& CDP_file_postfix);

        bool showProgressBar() const;
        void setShowProgressBar(bool show_progress_bar);

        const std::string& tag() const;
        void setTag(const std::string& tag);

        Scaling scaling() const;
        void setScaling(const Scaling& scaling);

        std::string workingDir() const;
        void setWorkingDir(const std::string& working_dir);

    private:
        bool m_show_progress_bar;
        bool m_gzip_cdp;
        bool m_gzip_mass;
        int m_mass_log_interval;
        int m_cdp_log_interval;
        std::string m_mass_file_postfix;
        std::string m_cdp_file_postfix;
        std::string m_tag;
        Scaling m_scaling;
        std::string m_working_dir;
    };

    std::string toString(LogParameter::Scaling scaling);
    LogParameter::Scaling scalingFromString(const std::string& str, bool* ok = nullptr);

    // Overall paramter pack
    class Parameter
    {
      public:
        Parameter();
        bool isValid(std::string& error_string) const;
        std::string overviewString() const;

        SystemParameter& systemParameter();
        const SystemParameter& systemParameter() const;

        PKParameter& pkParameter();
        const PKParameter& pkParameter() const;

        const SinkParameter& sinkParameter() const;
        SinkParameter& sinkParameter();

        const VehicleParameter& vehicleParameter() const;
        VehicleParameter& vehicleParameter();

        const std::vector<LayerParameter>& layer_parameter() const;
        const LayerParameter& layer(int idx) const;
        LayerParameter& layer(int idx);
        int layerCount() const;
        void addLayer(const LayerParameter& parameter);

        const LogParameter& logParameter() const;
        LogParameter& logParameter();

    private:
        SystemParameter m_system_parameter;
        LogParameter m_log_parameter;
        PKParameter m_pk_parameter;
        SinkParameter m_sink_parameter;
        VehicleParameter m_vehicle_parameter;
        std::vector<LayerParameter> m_layer_parameter;
    };
}

#endif  // PARAMETER_H
