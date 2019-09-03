#include "versioninfo.h"
#include <sstream>
#include <ctime>

namespace sc
{
    VersionInfo::VersionInfo() : m_maj_ver(0), m_min_ver(0), m_patch_level(0)
    {
    }

    VersionInfo::VersionInfo(const std::string& app_name, int maj_ver, int min_ver, int patch_lvl)
        : m_appname(app_name), m_maj_ver(maj_ver), m_min_ver(min_ver), m_patch_level(patch_lvl)
    {
    }

    const std::string& VersionInfo::appName() const
    {
        return m_appname;
    }

    void VersionInfo::setAppName(const std::string& app_name)
    {
        m_appname = app_name;
    }

    std::string VersionInfo::fullName() const
    {
        return m_full_name;
    }

    void VersionInfo::setFullName(const std::string& full_name)
    {
        m_full_name = full_name;
    }

    int VersionInfo::majVer() const
    {
        return m_maj_ver;
    }

    void VersionInfo::setMajVer(int maj_ver)
    {
        m_maj_ver = maj_ver;
    }

    int VersionInfo::minVer() const
    {
        return m_min_ver;
    }

    void VersionInfo::setMinVer(int min_ver)
    {
        m_min_ver = min_ver;
    }

    int VersionInfo::patchLevel() const
    {
        return m_patch_level;
    }

    void VersionInfo::setPatchLevel(int patch_level)
    {
        m_patch_level = patch_level;
    }

    std::string VersionInfo::buildId() const
    {
        return m_build_id;
    }

    void VersionInfo::setBuildId(const std::string& value)
    {
        m_build_id = value;
    }

    std::string VersionInfo::versionString() const
    {
        std::stringstream ss;

        ss << this->majVer() << "." << this->minVer() << "." << this->patchLevel();
        ss << " (Build id: " << this->buildId() << ")\n";

        return ss.str();
    }

    std::string VersionInfo::copyrightNote() const
    {
        return m_copyright_note;
    }

    void VersionInfo::setCopyrightNote(const std::string& copyright_note)
    {
        m_copyright_note = copyright_note;
    }

    std::ostream& operator<<(std::ostream& out, const VersionInfo& info)
    {
        const auto t_now = std::time(0);
        auto t_st = std::localtime(&t_now);
        const auto current_year = t_st->tm_year + 1900;

        out << info.appName() << " - " << info.fullName() << "\n";
        out << "(c) " << info.copyrightNote() << " (" << current_year << ")\n";
        out << "Version  : " << info.majVer() << "." << info.minVer() << "." << info.patchLevel()
            << "\n";
        out << "Build id : " << info.buildId() << "\n";

        return out;
    }
}
