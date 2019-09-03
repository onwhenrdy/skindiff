#include "cmdlineparser.h"

#include "jsonparser.h"
#include <fstream>
#include <iostream>
#include <sstream>

namespace sc
{
    CmdLineParser::CmdLineParser()
    {
    }

    CmdLineParser::Status CmdLineParser::parse(const std::vector<std::string>& parameter)
    {
        // clear
        m_parameter         = Parameter();
        const auto n_params = parameter.size();

        // check for template and/or config files
        if (n_params == 1)
        {
            if (parameter[0] == "--template")
            {
                return CmdLineParser::Status::Write_CFG_Template;
            }
            else if (parameter[0] == "--version")
            {
                return CmdLineParser::Status::Version_Info;
            }
            else
            {
                // parse config file
                JsonParser parser;
                const auto ok = parser.parseFromFile(parameter[0]);
                if (!ok)
                {
                    m_last_error = parser.lastError();
                    return CmdLineParser::Status::Error;
                }
                else
                {
                    m_parameter = parser.parameter();
                    return CmdLineParser::Status::Parsed_CFG_File;
                }
            }
        }

        // the minimum number is 19 parameter
        // the last parameter is always the string tag

        if (n_params < 19)
        {
            m_last_error = "Need at least 19 input parameters.";
            return CmdLineParser::Status::Error;
        }
        if (n_params != 19 && n_params != 20 && n_params != 21 && n_params != 23)
        {
            m_last_error = "Need 19, 20, 21 or 23 input parameters.";
            return CmdLineParser::Status::Error;
        }

        // start parsing
        LayerParameter sc;
        sc.setName("SC");
        LayerParameter dsl;
        dsl.setName("DSL");

        bool ok           = false;
        const auto c_init = this->strToDouble(parameter[0], &ok);
        if (!ok)
        {
            m_last_error = "C_0 is not a double.";
            return CmdLineParser::Status::Error;
        }
        if (c_init < 0.0)
        {
            m_last_error = "C_0 < 0.0";
            return CmdLineParser::Status::Error;
        }
        m_parameter.vehicleParameter().setCInit(c_init);

        const auto d_donor = this->strToDouble(parameter[1], &ok);
        if (!ok)
        {
            m_last_error = "D_Donor is not a double.";
            return CmdLineParser::Status::Error;
        }
        if (d_donor <= 0.0)
        {
            m_last_error = "D_Donor <= 0.0";
            return CmdLineParser::Status::Error;
        }
        m_parameter.vehicleParameter().setD(d_donor);

        const auto d_sc = this->strToDouble(parameter[2], &ok);
        if (!ok)
        {
            m_last_error = "D_SC is not a double.";
            return CmdLineParser::Status::Error;
        }
        if (d_sc <= 0.0)
        {
            m_last_error = "D_SC <= 0.0";
            return CmdLineParser::Status::Error;
        }
        sc.setD(d_sc);

        const auto d_dsl = this->strToDouble(parameter[3], &ok);
        if (!ok)
        {
            m_last_error = "D_DSL is not a double.";
            return CmdLineParser::Status::Error;
        }
        if (d_dsl <= 0.0)
        {
            m_last_error = "D_DSL <= 0.0";
            return CmdLineParser::Status::Error;
        }
        dsl.setD(d_dsl);

        const auto k_sc = this->strToDouble(parameter[4], &ok);
        if (!ok)
        {
            m_last_error = "K_SC/Don is not a double.";
            return CmdLineParser::Status::Error;
        }
        if (k_sc <= 0.0)
        {
            m_last_error = "K_SC/Don <= 0.0";
            return CmdLineParser::Status::Error;
        }
        sc.setK(k_sc);

        const auto k_acc = this->strToDouble(parameter[5], &ok);
        if (!ok)
        {
            m_last_error = "K_DSL/Don is not a double.";
            return CmdLineParser::Status::Error;
        }
        if (k_acc <= 0.0)
        {
            m_last_error = "K_DSL/Don <= 0.0";
            return CmdLineParser::Status::Error;
        }
        dsl.setK(k_acc);

        const auto A = this->strToDouble(parameter[6], &ok);
        if (!ok)
        {
            m_last_error = "App Area is not a double.";
            return CmdLineParser::Status::Error;
        }
        if (A <= 0.0)
        {
            m_last_error = "App Area < 0.0";
            return CmdLineParser::Status::Error;
        }
        m_parameter.vehicleParameter().setAppArea(A);

        const auto l_cross_section = this->strToDouble(parameter[7], &ok);
        if (!ok)
        {
            m_last_error = "Lipid CS is not a double.";
            return CmdLineParser::Status::Error;
        }
        if (l_cross_section <= 0.0 || l_cross_section > 1.0)
        {
            m_last_error = "Lipid CS <= 0.0 or > 1.0";
            return CmdLineParser::Status::Error;
        }
        sc.setCrossSection(l_cross_section);

        const auto dsl_cross_section = this->strToDouble(parameter[8], &ok);
        if (!ok)
        {
            m_last_error = "DSL CS is not a double.";
            return CmdLineParser::Status::Error;
        }
        if (dsl_cross_section <= 0.0 || dsl_cross_section > 1.0)
        {
            m_last_error = "DSL CS <= 0.0 or > 1.0";
            return CmdLineParser::Status::Error;
        }
        dsl.setCrossSection(dsl_cross_section);

        const auto h_don = this->strToInt(parameter[9], &ok);
        if (!ok)
        {
            m_last_error = "h_Donor is not an int.";
            return CmdLineParser::Status::Error;
        }
        if (h_don < 1)
        {
            m_last_error = "h_Donor < 1";
            return CmdLineParser::Status::Error;
        }
        m_parameter.vehicleParameter().setHeight(h_don);

        const auto h_sc = this->strToInt(parameter[10], &ok);
        if (!ok)
        {
            m_last_error = "h_SC is not an int.";
            return CmdLineParser::Status::Error;
        }
        if (h_sc < 1)
        {
            m_last_error = "h_SC < 1";
            return CmdLineParser::Status::Error;
        }
        sc.setHeight(h_sc);

        const auto h_dsl = this->strToInt(parameter[11], &ok);
        if (!ok)
        {
            m_last_error = "h_DSL is not an int.";
            return CmdLineParser::Status::Error;
        }
        if (h_dsl < 1)
        {
            m_last_error = "h_DSL < 1";
            return CmdLineParser::Status::Error;
        }
        dsl.setHeight(h_dsl);

        const auto sim_time = this->strToInt(parameter[12], &ok);
        if (!ok)
        {
            m_last_error = "Sim time is not an int.";
            return CmdLineParser::Status::Error;
        }
        if (sim_time < 1)
        {
            m_last_error = "Sim time < 1";
            return CmdLineParser::Status::Error;
        }
        m_parameter.systemParameter().setSimulationTime(sim_time);

        const auto resolution = this->strToInt(parameter[13], &ok);
        if (!ok)
        {
            m_last_error = "Resolution is not an int.";
            return CmdLineParser::Status::Error;
        }
        if (resolution < 1)
        {
            m_last_error = "Resolution < 1";
            return CmdLineParser::Status::Error;
        }
        m_parameter.systemParameter().setResolution(resolution);

        const auto scaling_str = parameter[14];
        m_parameter.logParameter().setScaling(scalingFromString(scaling_str, &ok));
        if (!ok)
        {
            m_last_error = "Unknown scaling string: '" + scaling_str + "'";
            return CmdLineParser::Status::Error;
        }

        const auto disc_str = parameter[15];
        m_parameter.systemParameter().setDiscMethod(discMethodFromString(disc_str, &ok));
        if (!ok)
        {
            m_last_error = "Unknown discretization method: '" + disc_str + "'";
            return CmdLineParser::Status::Error;
        }

        const auto mb_str = parameter[16];
        m_parameter.systemParameter().setMatrixBuilderMethod(mbMethodFromString(mb_str, &ok));
        if (!ok)
        {
            m_last_error = "Unknown matrix builder method: '" + mb_str + "'";
            return CmdLineParser::Status::Error;
        }

        const auto dose_str = parameter[17];
        ok                  = (dose_str == "yes" || dose_str == "no");
        if (!ok)
        {
            m_last_error = "Unknown finite dose string (expected yes or no): '" + dose_str + "'";
            return CmdLineParser::Status::Error;
        }
        m_parameter.vehicleParameter().setFiniteDose((dose_str == "yes"));

        int i = 18;
        if (n_params >= 20)
        {
            const auto rem_at = this->strToInt(parameter[i]);
            if (!ok)
            {
                m_last_error = "Remove at is not an int.";
                return CmdLineParser::Status::Error;
            }
            if (rem_at < 0)
            {
                m_last_error = "Remove at < 0";
                return CmdLineParser::Status::Error;
            }
            m_parameter.vehicleParameter().setRemoveAt(rem_at);
            i++;
        }

        if (n_params >= 21)
        {
            const auto rep_after = this->strToInt(parameter[i]);
            if (!ok)
            {
                m_last_error = "Replicate after is not an int.";
                return CmdLineParser::Status::Error;
            }
            if (rep_after < 0)
            {
                m_last_error = "Replicate after < 0";
                return CmdLineParser::Status::Error;
            }

            m_parameter.vehicleParameter().setReplaceAfter(rep_after);
            i++;
        }

        if (n_params == 23)
        {
            m_parameter.pkParameter().setEnabled(true);
            const auto vd = this->strToDouble(parameter[i]);
            if (!ok)
            {
                m_last_error = "Vd is not a double.";
                return CmdLineParser::Status::Error;
            }
            if (vd <= 0)
            {
                m_last_error = "Vd <= 0";
                return CmdLineParser::Status::Error;
            }
            m_parameter.sinkParameter().setVd(vd);
            i++;

            const auto t_half = this->strToDouble(parameter[i]);
            if (!ok)
            {
                m_last_error = "t_half is not a double.";
                return CmdLineParser::Status::Error;
            }
            if (t_half <= 0.0)
            {
                m_last_error = "t_half <= 0";
                return CmdLineParser::Status::Error;
            }
            m_parameter.pkParameter().setThalf(t_half);
            i++;
        }

        const auto file_str = parameter[i];
        if (file_str.empty())
        {
            m_last_error = "File Tag is empty";
            return CmdLineParser::Status::Error;
        }
        m_parameter.logParameter().setTag(file_str);

        m_parameter.addLayer(sc);
        m_parameter.addLayer(dsl);

        return CmdLineParser::Status::Parsed_CMD_Line;
    }

    double CmdLineParser::strToDouble(const std::string& value, bool* ok)
    {
        double res = 0.0;
        if (ok)
        {
            *ok = true;
        }

        try
        {
            res = std::stod(value);
        }
        catch (...)
        {
            if (ok)
            {
                *ok = false;
            }
        }
        return res;
    }

    int CmdLineParser::strToInt(const std::string& value, bool* ok)
    {
        int res = 0;
        if (ok)
        {
            *ok = true;
        }

        try
        {
            res = std::stoi(value);
        }
        catch (...)
        {
            if (ok)
            {
                *ok = false;
            }
        }
        return res;
    }

    const Parameter& CmdLineParser::parameter() const
    {
        return m_parameter;
    }

    std::string CmdLineParser::cmdlineOptions()
    {
        std::stringstream ss;

        ss << "Option List:\n";
        ss << "--------------------------------\n";
        ss << "--template     : Creates a config file template\n";
        ss << "--version      : Outputs version information\n";
        ss << "[1] FILENAME   : Uses the config file parameters (ignores other parameter)\n\n";

        ss << "Cmdline Parameter List:\n";
        ss << "--------------------------------\n";
        ss << "[1]  C_0         [mg/ml]"
           << "\n";
        ss << "[2]  D_Donor     [um^2/min]"
           << "\n";
        ss << "[3]  D_SC        [um^2/min]"
           << "\n";
        ss << "[4]  D_DSL       [um^2/min]"
           << "\n";
        ss << "[5]  K_SC/Don    [no unit]"
           << "\n";
        ss << "[6]  K_DSL/Don   [no unit]"
           << "\n";
        ss << "[7]  App area    [cm^2]"
           << "\n";
        ss << "[8]  Lipid CS    ]0..1]"
           << "\n";
        ss << "[9]  DSL CS      ]0..1]"
           << "\n";
        ss << "[10] h_Donor     [um]"
           << "\n";
        ss << "[11] h_SC        [um]"
           << "\n";
        ss << "[12] h_DSL       [um]"
           << "\n";
        ss << "[13] Sim time    [min]"
           << "\n";
        ss << "[14] Resolution  [1/x um]"
           << "\n";
        ss << "[15] Scaling     [mg/ug/ng]"
           << "\n";
        ss << "[16] Disc method [EQUIDIST/BK]"
           << "\n";
        ss << "[17] MB method   [DSkin_1_3/DSkin_1_4]"
           << "\n";
        ss << "[18] Finite dose [yes/no]"
           << "\n";
        ss << "[19] Remove at   [min] (optional; 0 to disable)"
           << "\n";
        ss << "[20] Repl. after [min; interval] (optional; 0 to disable)"
           << "\n";
        ss << "[21] Vd          [ml]  (optional; 0 to disable - enables PK - t_half nedded!)"
           << "\n";
        ss << "[22] t_half      [min] (optional; 0 to disable - enables PK - Vd nedded!)"
           << "\n";
        ss << "[..] File tag    [string] (not optional!)"
           << "\n";

        return ss.str();
    }

    std::string CmdLineParser::lastError() const
    {
        return "Parsing error: " + m_last_error + "\n";
    }
}
