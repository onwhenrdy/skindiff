#ifndef SC_JSONPARSER_H
#define SC_JSONPARSER_H

#include <string>
#include "json.h"
#include "parameter.h"

namespace sc
{
    class JsonParser
    {
      public:
        JsonParser();
        static std::string templateString();
        bool parseFromFile(const std::string& filename);
        bool parseFromString(const std::string& parameter);
        std::string lastError() const;
        const Parameter& parameter() const;

    private:
        bool parseJson(const nlohmann::json& object);

    private:
        static const std::string m_template_string;
        std::string m_last_error;
        Parameter m_parameter;
    };
}

#endif  // JSONPARSER_H
