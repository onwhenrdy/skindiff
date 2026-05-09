#include "parameter.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <string>

namespace sc
{
    std::string_view toString(Scaling s) noexcept
    {
        switch (s)
        {
            case Scaling::MG: return "mg";
            case Scaling::UG: return "ug";
            case Scaling::NG: return "ng";
        }
        return "mg";
    }

    std::optional<Scaling> scalingFromString(std::string_view str) noexcept
    {
        std::string upper(str);
        std::transform(upper.begin(), upper.end(), upper.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        if (upper == "MG") return Scaling::MG;
        if (upper == "UG") return Scaling::UG;
        if (upper == "NG") return Scaling::NG;
        return std::nullopt;
    }

    double scaleFactor(Scaling s) noexcept
    {
        switch (s)
        {
            case Scaling::MG: return 1.0;
            case Scaling::UG: return 1.0e3;
            case Scaling::NG: return 1.0e6;
        }
        return 1.0;
    }

    namespace
    {
        std::optional<std::string> validate(const VehicleParams& v)
        {
            if (v.name.empty())   return "vehicle.name is empty";
            if (v.c_init  < 0.0)  return "vehicle.c_init < 0";
            if (v.app_area <= 0.0) return "vehicle.app_area <= 0";
            if (v.D       < 0.0)  return "vehicle.D < 0";
            if (v.height  < 3)    return "vehicle.height < 3 um";
            if (v.remove_at    < 0) return "vehicle.remove_at < 0";
            if (v.replace_after < 0) return "vehicle.replace_after < 0";
            return std::nullopt;
        }

        std::optional<std::string> validate(const LayerParams& l, std::size_t idx)
        {
            std::ostringstream tag;
            tag << "layer[" << idx << "].";
            if (l.name.empty())                          return tag.str() + "name is empty";
            if (l.c_init < 0.0)                          return tag.str() + "c_init < 0";
            if (l.D      < 0.0)                          return tag.str() + "D < 0";
            if (l.K      <= 0.0)                         return tag.str() + "K <= 0";
            if (l.cross_section <= 0.0 || l.cross_section > 1.0)
                return tag.str() + "cross_section not in (0, 1]";
            if (l.height < 3)                            return tag.str() + "height < 3 um";
            return std::nullopt;
        }

        std::optional<std::string> validate(const SinkParams& s)
        {
            if (s.name.empty())  return "sink.name is empty";
            if (s.Vd     <= 0.0) return "sink.Vd <= 0";
            if (s.c_init <  0.0) return "sink.c_init < 0";
            return std::nullopt;
        }

        std::optional<std::string> validate(const PKParams& pk)
        {
            if (pk.enabled && pk.thalf <= 0.0) return "pk.thalf <= 0";
            return std::nullopt;
        }

        std::optional<std::string> validate(const SystemParams& s)
        {
            if (s.resolution      <= 0)             return "sys.resolution <= 0";
            if (s.max_module      <= 0.0)           return "sys.max_module <= 0";
            if (s.simulation_time <= 0)             return "sys.simulation_time <= 0";
            return std::nullopt;
        }

        std::optional<std::string> validate(const LogParams& l)
        {
            if (l.mass_log_interval <= 0) return "log.mass_log_interval <= 0";
            if (l.cdp_log_interval  <= 0) return "log.cdp_log_interval <= 0";
            return std::nullopt;
        }
    }

    std::optional<std::string> validate(const Parameters& p)
    {
        if (auto err = validate(p.sys))     return err;
        if (auto err = validate(p.log))     return err;
        if (auto err = validate(p.pk))      return err;
        if (auto err = validate(p.sink))    return err;
        if (auto err = validate(p.vehicle)) return err;
        for (std::size_t i = 0; i < p.layers.size(); ++i)
        {
            if (auto err = validate(p.layers[i], i)) return err;
        }
        if (p.vehicle.removed() && p.layers.empty())
        {
            return "cannot remove the vehicle if no layers are defined";
        }
        return std::nullopt;
    }
}
