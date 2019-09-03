#include "jsonparser.h"
#include "geometry.h"
#include "matrixbuilder.h"
#include <fstream>
#include <iostream>

const std::string sc::JsonParser::m_template_string =
    R"({
    "sys" :
    {
        "disc_scheme" : "BK",
        "mb_method" : "DSkin_1_4",
        "resolution" : 1,
        "max_module" : 50.0,
        "mb_eta" : 0.6,
        "sim_time" : 600
    },

    "log" :
    {
        "file_tag" : "test",
        "mass_file_postfix" : "mass",
        "mass_file_gzip" : false,
        "cdp_file_postfix" : "cdp",
        "cdp_file_gzip" : true,
        "mass_log_interval" : 1,
        "cdp_log_interval" : 1,
        "scaling" : "mg",
        "show_progress" : true,
        "working_dir" : ""
    },

    "PK" :
    {
        "enabled" : true,
        "t_half" : 1.0
    },

    "compartments" :
    {
        "vehicle" :
        {
            "name" : "Donor",
            "finite_dose" : true,
            "c_init" : 1.0,
            "app_area" : 1.0,
            "h" : 30,
            "D" : 1.0,
            "replace_after" : 200,
            "remove_after" : 400,
            "log" : true,
            "log_cdp" : true
        },

        "sink" :
        {
            "name" : "Sink",
            "log" : true,
            "c_init" : 0.0,
            "Vd" : 1.0
        },

        "layers" :
        [
            {
                "name" : "SC",
                "log" : true,
                "log_cdp" : true,
                "c_init" : 0.0,
                "cross_section" : 1.0,
                "h" : 10,
                "D" : 1.0,
                "K" : 1.0
            },

            {
                "name" : "DSL",
                "log" : true,
                "log_cdp" : true,
                "c_init" : 0.0,
                "cross_section" : 1.0,
                "h" : 10,
                "D" : 1.0,
                "K" : 1.0
            }
        ]
    }
})";

using Json = nlohmann::json;

namespace sc
{
    JsonParser::JsonParser()
    {
    }

    std::string JsonParser::templateString()
    {
        return m_template_string;
    }

    bool JsonParser::parseFromFile(const std::string& filename)
    {
        std::ifstream file(filename);
        if (!file)
        {
            m_last_error = "Could not open file " + filename;
            return false;
        }

        Json j;
        try
        {
            file >> j;
            return this->parseJson(j);
        }
        catch (Json::parse_error& error)
        {
            m_last_error = std::string(error.what());
            return false;
        }

        return true;
    }

    bool JsonParser::parseFromString(const std::string& parameter)
    {
        Json j;
        try
        {
            j = Json::parse(parameter);
            return this->parseJson(j);
        }
        catch (Json::parse_error& error)
        {
            m_last_error = std::string(error.what());
            return false;
        }

        return true;
    }

    std::string JsonParser::lastError() const
    {
        return "Error parsing config string: " + m_last_error;
    }

    bool JsonParser::parseJson(const nlohmann::json& object)
    {
        m_parameter = Parameter();
        m_last_error.clear();

        // system
        const auto system = object.find("sys");
        if (system == object.end())
        {
            m_last_error = "Could not find <sys> section.";
            return false;
        }
        else
        {
            const auto sys_value = system.value();

            auto& params = m_parameter.systemParameter();
            params.setResolution(sys_value.value("resolution", 1));
            params.setSimulationTime(sys_value.value("sim_time", 60));
            params.setMaxModule(sys_value.value("max_module", 50.0));
            params.setEta(sys_value.value("mb_eta", 0.6));

            bool ok = false;
            params.setDiscMethod(
                discMethodFromString(sys_value.value("disc_scheme", "EQUIDIST"), &ok));
            if (!ok)
            {
                m_last_error = "Unknown disc_scheme found.";
                return false;
            }

            ok = false;
            params.setMatrixBuilderMethod(
                mbMethodFromString(sys_value.value("mb_method", "DSkin_1_5"), &ok));
            if (!ok)
            {
                m_last_error = "Unknown mb_method found.";
                return false;
            }
        }

        // log
        const auto log = object.find("log");
        if (log != object.end())
        {
            const auto log_value = log.value();
            auto& params         = m_parameter.logParameter();

            // all optional
            params.setTag(log_value.value("file_tag", "unknown"));
            params.setShowProgressBar(log_value.value("show_progress", true));

            bool ok = false;
            params.setScaling(scalingFromString(log_value.value("scaling", "mg"), &ok));
            if (!ok)
            {
                m_last_error = "Unknown scaling found.";
                return false;
            }

            params.setWorkingDir(log_value.value("working_dir", ""));
            params.setMassFilePostfix(log_value.value("mass_file_postfix", "mass"));
            params.setGzipMass(log_value.value("mass_file_gzip", false));
            params.setCDPFilePostfix(log_value.value("cdp_file_postfix", "cdp"));
            params.setGzipCDP(log_value.value("cdp_file_gzip", true));
            params.setMassLogInterval(log_value.value("mass_log_interval", 1));
            params.setCDPLogInterval(log_value.value("cdp_log_interval", 1));
        }

        // PK
        const auto pk = object.find("PK");
        if (pk != object.end())
        {
            const auto pk_value = pk.value();
            auto& params        = m_parameter.pkParameter();
            params.setEnabled(pk_value.value("enabled", true));
            if (pk_value.count("t_half") == 0)
            {
                m_last_error = "PK parameters need a t_half value.";
                return false;
            }
            params.setThalf(pk_value["t_half"]);
        }

        // compartments
        const auto comps = object.find("compartments");
        if (comps == object.end())
        {
            m_last_error = "Could not find <compartments> section.";
            return false;
        }
        const auto comps_value = comps.value();

        // sink (optional)
        const auto sink = comps_value.find("sink");
        if (sink != comps_value.end())
        {
            const auto sink_value = sink.value();
            m_parameter.sinkParameter().setName(sink_value.value("name", "Sink"));
            m_parameter.sinkParameter().setLog(sink_value.value("log", true));
            m_parameter.sinkParameter().setCInit(sink_value.value("c_init", 0.0));
            m_parameter.sinkParameter().setVd(sink_value.value("Vd", 1.0));
        }

        // vehicle (not optional)
        const auto vehicle = comps_value.find("vehicle");
        if (vehicle != comps_value.end())
        {
            const auto vehicle_value = vehicle.value();
            auto& params             = m_parameter.vehicleParameter();

            // optinals
            params.setAppArea(vehicle_value.value("app_area", 1.0));
            params.setName(vehicle_value.value("name", "Vehicle"));
            params.setLog(vehicle_value.value("log", true));
            params.setLogCDP(vehicle_value.value("log_cdp", false));
            params.setReplaceAfter(vehicle_value.value("replace_after", 0));
            params.setRemoveAt(vehicle_value.value("remove_after", 0));
            params.setFiniteDose(vehicle_value.value("finite_dose", true));

            // need at least c_init, h and D
            if (vehicle_value.count("c_init") == 0 || vehicle_value.count("h") == 0 ||
                vehicle_value.count("D") == 0)
            {
                m_last_error = "Vehicle section needs at least values for c_init, h and D.";
                return false;
            }
            params.setCInit(vehicle_value["c_init"]);
            params.setD(vehicle_value["D"]);
            params.setHeight(vehicle_value["h"]);
        }

        // layers
        const auto layers = comps_value.find("layers");
        if (layers != comps_value.end())
        {
            const auto& layers_value = layers.value();
            if (!layers_value.is_array())
            {
                m_last_error = "Layers definition is malformated. Expected an array.";
                return false;
            }

            for (const auto& l : layers_value)
            {
                LayerParameter param;
                // optionals
                param.setLog(l.value("log", true));
                param.setLogCDP(l.value("log_cdp", false));
                param.setCrossSection(l.value("cross_section", 1.0));
                param.setCInit(l.value("c_init", 0.0));

                // need at least a name, h, D and K
                if (l.count("name") == 0 || l.count("h") == 0 || l.count("D") == 0 ||
                    l.count("K") == 0)
                {
                    m_last_error = "Layers need at least values for name, h, D and K.";
                    return false;
                }

                param.setD(l["D"]);
                param.setK(l["K"]);
                param.setHeight(l["h"]);
                param.setName(l["name"]);

                m_parameter.addLayer(param);
            }
        }

        const auto ok = m_parameter.isValid(m_last_error);
        return ok;
    }

    const Parameter& JsonParser::parameter() const
    {
        return m_parameter;
    }
}
