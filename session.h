#ifndef SC_SESSION_H
#define SC_SESSION_H

#include <string>
#include "versioninfo.h"

namespace sc
{
    class Session
    {
      public:
        Session(const VersionInfo& v_info, const std::string& id ="no_id");

        std::string id() const;
        void setId(const std::string &id);

        bool runFromCmdLineArgs(int nargs, const char* args[]);

        bool warningsShown() const;
        void showWarnings(bool show_warnings);

        bool infosShown() const;
        void showInfos(bool show_infos);

        const VersionInfo& versionInfo() const;

    private:
        std::string m_id;
        bool m_show_infos;
        bool m_show_warnings;
        VersionInfo m_version_info;
    };
}
#endif  // SESSION_H
