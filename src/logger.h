#ifndef SC_LOGGER_H
#define SC_LOGGER_H

#include <cmath>
#include <vector>

namespace sc
{
    // Time-series of a scalar (mass per compartment, mass per sink).
    struct MassSeries
    {
        bool   enabled      = true;
        int    log_interval = 1;            // minutes
        std::vector<double> times;          // minutes
        std::vector<double> values;         // in user-selected scaling unit

        void reserve_for_total(int total_minutes)
        {
            const auto cap =
                1u + static_cast<unsigned>(std::floor(total_minutes / std::max(1, log_interval)));
            times.reserve(cap);
            values.reserve(cap);
        }

        void record(double t, double value)
        {
            times.push_back(t);
            values.push_back(value);
        }

        [[nodiscard]] bool should_log(double t) const noexcept
        {
            if (!enabled) return false;
            return static_cast<int>(t) % log_interval == 0;
        }
    };

    // Time-series of a 1-D concentration profile (CDP) for a single compartment.
    // depths_um is fixed for the lifetime of the series; conc_per_time has one
    // entry per logged time-point, each of length depths_um.size().
    struct CdpSeries
    {
        bool   enabled      = false;
        int    log_interval = 1;            // minutes
        std::vector<double> depths_um;
        std::vector<double> times;          // minutes
        std::vector<std::vector<double>> conc_per_time;  // [time_idx][depth_idx]

        void reserve_for_total(int total_minutes)
        {
            const auto cap =
                1u + static_cast<unsigned>(std::floor(total_minutes / std::max(1, log_interval)));
            times.reserve(cap);
            conc_per_time.reserve(cap);
        }

        void record(double t, std::vector<double> profile)
        {
            times.push_back(t);
            conc_per_time.push_back(std::move(profile));
        }

        [[nodiscard]] bool should_log(double t) const noexcept
        {
            if (!enabled) return false;
            return static_cast<int>(t) % log_interval == 0;
        }
    };
}

#endif  // SC_LOGGER_H
