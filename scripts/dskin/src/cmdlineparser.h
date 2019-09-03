#ifndef SC_CMDLINEPARSER_H
#define SC_CMDLINEPARSER_H

#include "parameter.h"
#include <string>
#include <vector>

namespace sc
{
    class CmdLineParser
    {
      public:
        enum class Status
        {
            Error,
            Write_CFG_Template,
            Version_Info,
            Parsed_CFG_File,
            Parsed_CMD_Line
        };

        CmdLineParser();
        Status parse(const std::vector<std::string>& parameter);
        std::string lastError() const;

        const Parameter& parameter() const;

        static std::string cmdlineOptions();

      private:
        double strToDouble(const std::string& value, bool* ok = nullptr);
        int strToInt(const std::string& value, bool* ok = nullptr);

      private:
        std::string m_last_error;
        Parameter m_parameter;
    };
}
#endif  // CMDLINEPARSER_H
