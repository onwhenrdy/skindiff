#ifndef SC_PARAMETER_H
#define SC_PARAMETER_H

#include "geometry.h"
#include "matrixbuilder.h"

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace sc
{
    enum class Scaling
    {
        MG,
        UG,
        NG
    };

    [[nodiscard]] std::string_view toString(Scaling s) noexcept;
    [[nodiscard]] std::optional<Scaling> scalingFromString(std::string_view str) noexcept;
    [[nodiscard]] double scaleFactor(Scaling s) noexcept;

    struct VehicleParams
    {
        std::string name   = "Vehicle";
        double c_init      = 1.0;   // mg/ml
        double app_area    = 1.0;   // cm^2
        double D           = 1.0;   // um^2/min
        int    height      = 10;    // um
        int    replace_after = 0;   // min, 0 = disabled
        int    remove_at     = 0;   // min, 0 = disabled
        bool   finite_dose   = true;
        bool   log_mass      = true;
        bool   log_cdp       = false;

        [[nodiscard]] bool replaces() const noexcept { return replace_after > 0; }
        [[nodiscard]] bool removed() const noexcept { return remove_at > 0; }
    };

    struct LayerParams
    {
        std::string name;
        double c_init        = 0.0;   // mg/ml
        double D             = 1.0;   // um^2/min
        double K             = 1.0;   // partition coefficient relative to vehicle
        double cross_section = 1.0;   // (0, 1]
        int    height        = 10;    // um
        bool   log_mass      = true;
        bool   log_cdp       = false;
    };

    struct SinkParams
    {
        std::string name = "Sink";
        double c_init    = 0.0;   // mg/ml
        double Vd        = 1.0;   // ml
        bool   log_mass  = true;
    };

    struct SystemParams
    {
        int    resolution      = 1;     // sub-steps per um at the smallest-D compartment
        double max_module      = 50.0;  // sub-step stability target
        int    simulation_time = 600;   // min
    };

    struct LogParams
    {
        Scaling scaling           = Scaling::MG;
        int     mass_log_interval = 1;   // min
        int     cdp_log_interval  = 1;   // min
    };

    struct Parameters
    {
        SystemParams              sys;
        LogParams                 log;
        SinkParams                sink;
        VehicleParams             vehicle;
        std::vector<LayerParams>  layers;
    };

    // Returns std::nullopt on success, error message otherwise.
    [[nodiscard]] std::optional<std::string> validate(const Parameters& p);
}

#endif  // SC_PARAMETER_H
