#include "session.h"
#include "cmdlineparser.h"
#include "helper.h"
#include "jsonparser.h"
#include "systemcmd.h"
#include <fstream>
#include <iostream>

namespace sc
{
    Session::Session(const VersionInfo& v_info, const std::string& id)
        : m_id(id), m_show_infos(true), m_show_warnings(true), m_version_info(v_info)
    {
    }

    std::string Session::id() const
    {
        return m_id;
    }

    void Session::setId(const std::string& id)
    {
        m_id = id;
    }

    bool Session::runFromCmdLineArgs(int nargs, const char* args[])
    {
        // parse
        if (nargs < 2)
        {
            std::cout << CmdLineParser::cmdlineOptions();
            return true;
        }

        // parse cmline args
        std::vector<std::string> params(args + 1, args + nargs);
        CmdLineParser cmdparser;
        auto status = cmdparser.parse(params);
        if (status == CmdLineParser::Status::Error)
        {
            std::cerr << cmdparser.lastError();
            return false;
        }
        else if (status == CmdLineParser::Status::Version_Info)
        {
            std::cout << m_version_info.versionString();
            return true;
        }
        else if (status == CmdLineParser::Status::Write_CFG_Template)
        {
            static std::string cfg_filename("dskin_config.json");
            std::ofstream file(cfg_filename);
            if (file)
            {
                file << JsonParser::templateString();
                file.close();
                std::cout << "Wrote DSkin config template to file: " << cfg_filename << "\n";
            }
            else
            {
                std::cerr << "Could not write to file: " << cfg_filename << "\n";
                return false;
            }

            return true;
        }

        // print input parameter overview
        const auto& parameter = cmdparser.parameter();
        if (m_show_infos)
        {
            std::cout << m_version_info << "\n";
            std::cout << parameter.overviewString() << "\n";
        }

        SystemCmd system(parameter);
        auto result = system.run();
        auto cmdline_result = 1;
        if (result == System::Result::Executed || result == System::Result::Stopped)
        {
            auto ok = system.writeLogsToFiles();
            cmdline_result = (ok) ? 0 : 1;
        }

        if (m_show_infos)
        {
            std::cout << "\nComputation done.\n";
        }

        return cmdline_result;
    }

    bool Session::warningsShown() const
    {
        return m_show_warnings;
    }

    void Session::showWarnings(bool show_warnings)
    {
        m_show_warnings = show_warnings;
    }

    bool Session::infosShown() const
    {
        return m_show_infos;
    }

    void Session::showInfos(bool show_infos)
    {
        m_show_infos = show_infos;
    }

    const VersionInfo& Session::versionInfo() const
    {
        return m_version_info;
    }
}
