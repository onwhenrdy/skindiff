#include "parameter.h"
#include <algorithm>
#include <cctype>
#include <sstream>

namespace sc
{
    SystemParameter::SystemParameter()
        : m_disc_method(Geometry::DiscMethod::EQUI_DIST)
        , m_matrix_builder_method(MatrixBuilder::Method::DSkin_1_3)
        , m_resolution(1)
        , m_max_module(50.0)
        , m_eta(0.6)
        , m_simulation_time(60)
    {
    }

    bool SystemParameter::isValid(std::string& error_string) const
    {
        if (m_resolution <= 0)
        {
            error_string = "Resulution is <= 0.";
            return false;
        }

        if (m_max_module <= 0.0)
        {
            error_string = "Max module is <= 0.";
            return false;
        }

        if (m_simulation_time <= 0)
        {
            error_string = "Simulation time is <= 0.";
            return false;
        }

        if (m_eta <= 0 || m_eta > 1.0)
        {
            error_string = "mb_eta is <= 0 or > 1.0.";
            return false;
        }

        return true;
    }

    std::string SystemParameter::overviewString() const
    {
        std::stringstream ss;
        ss << "System Parameter:\n";
        ss << "--------------------------------\n";
        ss << "Discretization method  : " << toString(m_disc_method) << "\n";
        ss << "Matrix builder method  : " << toString(m_matrix_builder_method) << "\n\n";

        ss << "Sim time     [min]     : " << m_simulation_time << "\n";
        ss << "Resolution   [1/x um]  : " << m_resolution << "\n";
        ss << "MB scal. factor (eta)  : " << m_eta << "\n";
        ss << "Max Module             : " << m_max_module << "\n";

        return ss.str();
    }

    Geometry::DiscMethod SystemParameter::discMethod() const
    {
        return m_disc_method;
    }

    void SystemParameter::setDiscMethod(const Geometry::DiscMethod& disc_method)
    {
        m_disc_method = disc_method;
    }

    int SystemParameter::resolution() const
    {
        return m_resolution;
    }

    void SystemParameter::setResolution(int resolution)
    {
        m_resolution = resolution;
    }

    int SystemParameter::simulationTime() const
    {
        return m_simulation_time;
    }

    void SystemParameter::setSimulationTime(int simulation_time)
    {
        m_simulation_time = simulation_time;
    }

    double SystemParameter::maxModule() const
    {
        return m_max_module;
    }

    void SystemParameter::setMaxModule(double max_module)
    {
        m_max_module = max_module;
    }

    MatrixBuilder::Method SystemParameter::matrixBuilderMethod() const
    {
        return m_matrix_builder_method;
    }

    void SystemParameter::setMatrixBuilderMethod(const MatrixBuilder::Method& matrix_builder_method)
    {
        m_matrix_builder_method = matrix_builder_method;
    }

    double SystemParameter::eta() const
    {
        return m_eta;
    }

    void SystemParameter::setEta(double eta)
    {
        m_eta = eta;
    }

    Parameter::Parameter()
    {
    }

    bool Parameter::isValid(std::string& error_string) const
    {
        auto ok = m_system_parameter.isValid(error_string);
        if (!ok)
        {
            return false;
        }

        ok = m_log_parameter.isValid(error_string);
        if (!ok)
        {
            return false;
        }

        ok = m_pk_parameter.isValid(error_string);
        if (!ok)
        {
            return false;
        }

        ok = m_sink_parameter.isValid(error_string);
        if (!ok)
        {
            return false;
        }

        ok = m_vehicle_parameter.isValid(error_string);
        if (!ok)
        {
            return false;
        }

        for (const auto& layer : m_layer_parameter)
        {
            ok = layer.isValid(error_string);
            if (!ok)
            {
                return false;
            }
        }

        if (m_vehicle_parameter.remove() && this->layerCount() < 1)
        {
            error_string = "Cannot remove the vehicle if no layer is defined.";
            return false;
        }

        return ok;
    }

    std::string Parameter::overviewString() const
    {
        std::string result = m_system_parameter.overviewString();
        result += "\n";
        result += m_log_parameter.overviewString();
        result += "\n";
        result += m_pk_parameter.overviewString();
        result += "\n";
        result += m_vehicle_parameter.overviewString();
        result += "\n";
        result += m_sink_parameter.overviewString();
        for (const auto& layer : m_layer_parameter)
        {
            result += "\n";
            result += layer.overviewString();
        }

        return result;
    }

    SystemParameter& Parameter::systemParameter()
    {
        return m_system_parameter;
    }

    const SystemParameter& Parameter::systemParameter() const
    {
        return m_system_parameter;
    }

    PKParameter& Parameter::pkParameter()
    {
        return m_pk_parameter;
    }

    const PKParameter& Parameter::pkParameter() const
    {
        return m_pk_parameter;
    }

    const SinkParameter& Parameter::sinkParameter() const
    {
        return m_sink_parameter;
    }

    SinkParameter& Parameter::sinkParameter()
    {
        return m_sink_parameter;
    }

    VehicleParameter& Parameter::vehicleParameter()
    {
        return m_vehicle_parameter;
    }

    const std::vector<LayerParameter>& Parameter::layer_parameter() const
    {
        return m_layer_parameter;
    }

    const LayerParameter& Parameter::layer(int idx) const
    {
        return m_layer_parameter[idx];
    }

    LayerParameter& Parameter::layer(int idx)
    {
        return m_layer_parameter[idx];
    }

    int Parameter::layerCount() const
    {
        return m_layer_parameter.size();
    }

    void Parameter::addLayer(const LayerParameter& parameter)
    {
        m_layer_parameter.push_back(parameter);
    }

    LogParameter& Parameter::logParameter()
    {
        return m_log_parameter;
    }

    const LogParameter& Parameter::logParameter() const
    {
        return m_log_parameter;
    }

    const VehicleParameter& Parameter::vehicleParameter() const
    {
        return m_vehicle_parameter;
    }

    std::string toString(LogParameter::Scaling scaling)
    {
        switch (scaling)
        {
            case LogParameter::Scaling::MG:
                return "mg";
                break;
            case LogParameter::Scaling::UG:
                return "ug";
                break;
            case LogParameter::Scaling::NG:
                return "ng";
                break;
            default:
                return "unknown";
                break;
        }
    }

    LogParameter::Scaling scalingFromString(const std::string& str, bool* ok)
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
        if (upp_str == "MG")
        {
            return LogParameter::Scaling::MG;
        }

        if (upp_str == "NG")
        {
            return LogParameter::Scaling::NG;
        }

        if (upp_str == "UG")
        {
            return LogParameter::Scaling::UG;
        }

        if (ok)
        {
            *ok = false;
        }

        return LogParameter::Scaling::MG;
    }

    PKParameter::PKParameter() : m_enabled(false), m_thalf(0.0)
    {
    }

    bool PKParameter::isValid(std::string& error_string) const
    {
        if (m_enabled)
        {
            if (m_thalf <= 0.0)
            {
                error_string = "t_half <= 0.";
                return false;
            }
        }

        return true;
    }

    std::string PKParameter::overviewString() const
    {
        std::stringstream ss;
        ss << "PK Parameter:\n";
        ss << "--------------------------------\n";
        ss << "Enabled                : " << (m_enabled ? "yes" : "no") << "\n";
        if (m_enabled)
        {
            ss << "t 1/2       [h]        : " << m_thalf << "\n";
        }

        return ss.str();
    }

    bool PKParameter::enabled() const
    {
        return m_enabled;
    }

    void PKParameter::setEnabled(bool enabled)
    {
        m_enabled = enabled;
    }

    double PKParameter::thalf() const
    {
        return m_thalf;
    }

    void PKParameter::setThalf(double thalf)
    {
        m_thalf = thalf;
    }

    SinkParameter::SinkParameter() : m_log(true), m_name("Sink"), m_Vd(1.0), m_c_init(0.0)
    {
    }

    bool SinkParameter::isValid(std::string& error_string) const
    {
        if (m_name.empty())
        {
            error_string = "Sink name is empty.";
            return false;
        }

        if (m_Vd <= 0.0)
        {
            error_string = "Vd <= 0.0";
            return false;
        }

        if (m_c_init < 0.0)
        {
            error_string = "Sink C_init < 0.0";
            return false;
        }

        return true;
    }

    std::string SinkParameter::overviewString() const
    {
        std::stringstream ss;
        ss << "Sink Parameter:\n";
        ss << "--------------------------------\n";
        ss << "Name                   : " << m_name << "\n";
        ss << "Vd          [ml]       : " << m_Vd << "\n";
        ss << "C init      [mg/ml]    : " << m_c_init << "\n";
        ss << "Log Compartment        : " << (m_log ? "yes" : "no") << "\n";

        return ss.str();
    }

    bool SinkParameter::log() const
    {
        return m_log;
    }

    void SinkParameter::setLog(bool log)
    {
        m_log = log;
    }

    const std::string& SinkParameter::name() const
    {
        return m_name;
    }

    void SinkParameter::setName(const std::string& name)
    {
        m_name = name;
    }

    double SinkParameter::Vd() const
    {
        return m_Vd;
    }

    void SinkParameter::setVd(double Vd)
    {
        m_Vd = Vd;
    }

    double SinkParameter::cInit() const
    {
        return m_c_init;
    }

    void SinkParameter::setCInit(double value)
    {
        m_c_init = value;
    }

    VehicleParameter::VehicleParameter()
        : m_log(true)
        , m_log_cdp(false)
        , m_name("Vehicle")
        , m_c_init(1.0)
        , m_app_area(1.0)
        , m_D(1.0)
        , m_height(10)
        , m_replace_after(0)
        , m_remove_at(0)
        , m_finite_dose(true)
    {
    }

    bool VehicleParameter::isValid(std::string& error_string) const
    {
        if (m_name.empty())
        {
            error_string = "Vehicle name is empty.";
            return false;
        }

        if (m_c_init < 0.0)
        {
            error_string = "Vehicle C_init < 0.0.";
            return false;
        }

        if (m_app_area <= 0.0)
        {
            error_string = "Vehicle App Area <= 0.0.";
            return false;
        }

        if (m_D < 0.0)
        {
            error_string = "Vehicle D < 0.0.";
            return false;
        }

        if (m_height <= 2)
        {
            error_string = "Vehicle height < 2.";
            return false;
        }

        if (m_remove_at < 0)
        {
            error_string = "Vehicle remove at < 0.";
            return false;
        }

        if (m_replace_after < 0)
        {
            error_string = "Vehicle replace after < 0.";
            return false;
        }

        return true;
    }

    std::string VehicleParameter::overviewString() const
    {
        std::stringstream ss;
        ss << "Vehicle Parameter:\n";
        ss << "--------------------------------\n";
        ss << "Name                   : " << m_name << "\n";
        ss << "Log Mass               : " << (m_log ? "yes" : "no") << "\n";
        ss << "Log CDP                : " << (m_log_cdp ? "yes" : "no") << "\n";
        ss << "C init      [mg/ml]    : " << m_c_init << "\n";
        ss << "App Area    [cm^2]     : " << m_app_area << "\n";
        ss << "h           [um]       : " << m_height << "\n";
        ss << "D           [um^2/min] : " << m_D << "\n";
        ss << "Remove vehicle         : " << ((this->remove()) ? "yes" : "no") << "\n";
        if (this->remove())
        {
            ss << "Remove at   [min]      : " << m_remove_at << "\n";
        }

        ss << "Replace vehicle        : " << ((this->replace()) ? "yes" : "no") << "\n";
        if (this->replace())
        {
            ss << "Repl. after [min]      : " << m_replace_after << "\n";
        }
        ss << "Finite dose            : " << ((this->finiteDose()) ? "yes" : "no") << "\n";
        return ss.str();
    }

    bool VehicleParameter::log() const
    {
        return m_log;
    }

    void VehicleParameter::setLog(bool log)
    {
        m_log = log;
    }

    const std::string& VehicleParameter::name() const
    {
        return m_name;
    }

    void VehicleParameter::setName(const std::string& name)
    {
        m_name = name;
    }

    double VehicleParameter::cInit() const
    {
        return m_c_init;
    }

    double VehicleParameter::appArea() const
    {
        return m_app_area;
    }

    void VehicleParameter::setAppArea(double app_area)
    {
        m_app_area = app_area;
    }

    int VehicleParameter::height() const
    {
        return m_height;
    }

    void VehicleParameter::setHeight(int height)
    {
        m_height = height;
    }

    double VehicleParameter::D() const
    {
        return m_D;
    }

    void VehicleParameter::setD(double D)
    {
        m_D = D;
    }

    bool VehicleParameter::finiteDose() const
    {
        return m_finite_dose;
    }

    void VehicleParameter::setFiniteDose(bool finite_dose)
    {
        m_finite_dose = finite_dose;
    }

    bool VehicleParameter::logCDP() const
    {
        return m_log_cdp;
    }

    void VehicleParameter::setLogCDP(bool log_cdp)
    {
        m_log_cdp = log_cdp;
    }

    int VehicleParameter::replaceAfter() const
    {
        return m_replace_after;
    }

    void VehicleParameter::setReplaceAfter(int replace_after)
    {
        m_replace_after = replace_after;
    }

    bool VehicleParameter::replace() const
    {
        return m_replace_after > 0;
    }

    int VehicleParameter::removeAt() const
    {
        return m_remove_at;
    }

    void VehicleParameter::setRemoveAt(int remove_at)
    {
        m_remove_at = remove_at;
    }

    bool VehicleParameter::remove() const
    {
        return m_remove_at > 0;
    }

    void VehicleParameter::setCInit(double value)
    {
        m_c_init = value;
    }

    LayerParameter::LayerParameter()
        : m_log(true)
        , m_log_cdp(false)
        , m_name("")
        , m_c_init(0.0)
        , m_D(1.0)
        , m_K(1.0)
        , m_cross_section(1.0)
        , m_height(10)
    {
    }

    bool LayerParameter::isValid(std::string& error_string) const
    {
        if (m_name.empty())
        {
            error_string = "Layer name is empty.";
            return false;
        }

        if (m_c_init < 0.0)
        {
            error_string = "Layer C_init < 0.0.";
            return false;
        }

        if (m_D < 0.0)
        {
            error_string = "Layer D < 0.0.";
            return false;
        }

        if (m_K <= 0.0)
        {
            error_string = "Layer K <= 0.0.";
            return false;
        }

        if (m_cross_section <= 0.0 || m_cross_section > 1.0)
        {
            error_string = "Layer cross section not in ]0,1].";
            return false;
        }

        if (m_height <= 2)
        {
            error_string = "Layer height < 2.";
            return false;
        }

        return true;
    }

    std::string LayerParameter::overviewString() const
    {
        std::stringstream ss;
        ss << "Layer Parameter:\n";
        ss << "--------------------------------\n";
        ss << "Name                   : " << m_name << "\n";
        ss << "Log Mass               : " << (m_log ? "yes" : "no") << "\n";
        ss << "Log CDP                : " << (m_log_cdp ? "yes" : "no") << "\n";
        ss << "C init      [mg/ml]    : " << m_c_init << "\n";
        ss << "h           [um]       : " << m_height << "\n";
        ss << "D           [um^2/min] : " << m_D << "\n";
        ss << "K_Layer/Vehicle        : " << m_K << "\n";
        ss << "Layer CS    [%]        : " << m_cross_section * 100.0 << "\n";

        return ss.str();
    }

    bool LayerParameter::log() const
    {
        return m_log;
    }

    void LayerParameter::setLog(bool log)
    {
        m_log = log;
    }

    const std::string& LayerParameter::name() const
    {
        return m_name;
    }

    void LayerParameter::setName(const std::string& name)
    {
        m_name = name;
    }

    double LayerParameter::cInit() const
    {
        return m_c_init;
    }

    void LayerParameter::setCInit(double value)
    {
        m_c_init = value;
    }

    int LayerParameter::height() const
    {
        return m_height;
    }

    void LayerParameter::setHeight(int height)
    {
        m_height = height;
    }

    double LayerParameter::D() const
    {
        return m_D;
    }

    void LayerParameter::setD(double D)
    {
        m_D = D;
    }

    double LayerParameter::K() const
    {
        return m_K;
    }

    void LayerParameter::setK(double K)
    {
        m_K = K;
    }

    double LayerParameter::crossSection() const
    {
        return m_cross_section;
    }

    void LayerParameter::setCrossSection(double cross_section)
    {
        m_cross_section = cross_section;
    }

    bool LayerParameter::logCDP() const
    {
        return m_log_cdp;
    }

    void LayerParameter::setLogCDP(bool log_cdp)
    {
        m_log_cdp = log_cdp;
    }

    LogParameter::LogParameter()
        : m_show_progress_bar(true)
        , m_gzip_cdp(true)
        , m_gzip_mass(false)
        , m_mass_log_interval(1)
        , m_cdp_log_interval(1)
        , m_mass_file_postfix("mass")
        , m_cdp_file_postfix("cdp")
        , m_tag("unknonw")
        , m_scaling(LogParameter::Scaling::MG)
    {
    }

    bool LogParameter::isValid(std::string& error_string) const
    {
        if (m_mass_log_interval <= 0)
        {
            error_string = "Mass log interval <=0";
            return false;
        }

        if (m_cdp_log_interval <= 0)
        {
            error_string = "CDP log interval <=0";
            return false;
        }

        if (m_mass_file_postfix.empty())
        {
            error_string = "Mass file postfix is empty.";
            return false;
        }

        if (m_cdp_file_postfix.empty())
        {
            error_string = "CDP file postfix is empty.";
            return false;
        }

        if (m_tag.empty())
        {
            error_string = "File tag is empty.";
            return false;
        }

        return true;
    }

    std::string LogParameter::overviewString() const
    {
        std::stringstream ss;
        ss << "Log Parameter:\n";
        ss << "--------------------------------\n";
        ss << "File tag               : " << m_tag << "\n";
        ss << "Working directory      : " << m_working_dir << "\n";
        ss << "Mass logfile postfix   : " << m_mass_file_postfix << "\n";
        ss << "CDP logfile postfix    : " << m_cdp_file_postfix << "\n";
        ss << "Mass logfile gzip      : " << (m_gzip_mass ? "yes" : "no") << "\n";
        ss << "CDP logfile gzip       : " << (m_gzip_cdp ? "yes" : "no") << "\n";
        ss << "Mass log interv. [min] : " << m_mass_log_interval << "\n";
        ss << "CDP log interv. [min]  : " << m_cdp_log_interval << "\n";
        ss << "Scaling unit           : " << toString(m_scaling) << "\n";

        return ss.str();
    }

    bool LogParameter::gzipCDP() const
    {
        return m_gzip_cdp;
    }

    void LogParameter::setGzipCDP(bool gzip_CDP)
    {
        m_gzip_cdp = gzip_CDP;
    }

    bool LogParameter::gzipMass() const
    {
        return m_gzip_mass;
    }

    void LogParameter::setGzipMass(bool gzip_mass)
    {
        m_gzip_mass = gzip_mass;
    }

    int LogParameter::massLogInterval() const
    {
        return m_mass_log_interval;
    }

    void LogParameter::setMassLogInterval(int mass_log_interval)
    {
        m_mass_log_interval = mass_log_interval;
    }

    int LogParameter::CDPLogInterval() const
    {
        return m_cdp_log_interval;
    }

    void LogParameter::setCDPLogInterval(int cdp_log_interval)
    {
        m_cdp_log_interval = cdp_log_interval;
    }

    const std::string& LogParameter::massFilePostfix() const
    {
        return m_mass_file_postfix;
    }

    void LogParameter::setMassFilePostfix(const std::string& mass_file_postfix)
    {
        m_mass_file_postfix = mass_file_postfix;
    }

    const std::string& LogParameter::CDPFilePostfix() const
    {
        return m_cdp_file_postfix;
    }

    void LogParameter::setCDPFilePostfix(const std::string& CDP_file_postfix)
    {
        m_cdp_file_postfix = CDP_file_postfix;
    }

    bool LogParameter::showProgressBar() const
    {
        return m_show_progress_bar;
    }

    void LogParameter::setShowProgressBar(bool show_progress_bar)
    {
        m_show_progress_bar = show_progress_bar;
    }

    const std::string& LogParameter::tag() const
    {
        return m_tag;
    }

    void LogParameter::setTag(const std::string& tag)
    {
        m_tag = tag;
    }

    LogParameter::Scaling LogParameter::scaling() const
    {
        return m_scaling;
    }

    void LogParameter::setScaling(const Scaling& scaling)
    {
        m_scaling = scaling;
    }

    std::string LogParameter::workingDir() const
    {
        return m_working_dir;
    }

    void LogParameter::setWorkingDir(const std::string& working_dir)
    {
        m_working_dir = working_dir;
    }
}
