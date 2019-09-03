#ifndef SC_VERSIONINFO_H
#define SC_VERSIONINFO_H

#include <string>
#include <ostream>

namespace sc
{
    class VersionInfo
    {
      public:
        VersionInfo();
        VersionInfo(const std::string& app_name, int maj_ver, int min_ver, int patch_lvl);

        const std::string& appName() const;
        void setAppName(const std::string& app_name);

        std::string fullName() const;
        void setFullName(const std::string& fullName);

        int majVer() const;
        void setMajVer(int majVer);

        int minVer() const;
        void setMinVer(int minVer);

        int patchLevel() const;
        void setPatchLevel(int patchLevel);

        std::string buildId() const;
        void setBuildId(const std::string& value);

        std::string versionString() const;

        friend std::ostream& operator<<(std::ostream& out, const VersionInfo& info);

        std::string copyrightNote() const;
        void setCopyrightNote(const std::string& copyright_note);

    private:
        std::string m_appname;
        std::string m_full_name;
        std::string m_copyright_note;
        int m_maj_ver;
        int m_min_ver;
        int m_patch_level;
        std::string m_build_id;
    };

    std::ostream& operator<<(std::ostream& out, const VersionInfo& info);
}

#endif  // VERSIONINFO_H
