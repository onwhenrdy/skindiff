#ifndef SC_HELPER_H
#define SC_HELPER_H

#include <chrono>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace sc
{
    // Create string that can be used by R to plot/analyze the content
    template<typename T>
    std::string toRVector(const std::vector<T>& vals, const std::string& var_name = "a")
    {
        const int size = vals.size();

        std::stringstream ss;
        ss << var_name << " <- c(";
        if (size == 1)
        {
            ss << vals[0];
        }
        else if (size > 1)
        {
            for (int i = 0; i < size - 1; ++i)
            {
                ss << vals[i] << ", ";
            }
            ss << vals[size - 1];
        }
        ss << ")\n";

        return ss.str();
    }

    // small and simple benchmarking tool
    namespace bench
    {
        template<typename TimeT = std::chrono::milliseconds,
                 class ClockT   = std::chrono::system_clock>
        struct measure
        {
            template<typename F, typename... Args>
            static typename TimeT::rep execution(F&& func, Args&&... args)
            {
                auto start = ClockT::now();
                std::forward<decltype(func)>(func)(std::forward<Args>(args)...);
                auto duration = std::chrono::duration_cast<TimeT>(ClockT::now() - start);
                return duration.count();
            }

            template<typename F, typename... Args>
            static TimeT duration(F&& func, Args&&... args)
            {
                auto start = ClockT::now();
                std::forward<decltype(func)>(func)(std::forward<Args>(args)...);
                return std::chrono::duration_cast<TimeT>(ClockT::now() - start);
            }
        };
    }
}
#endif  // HELPER_H
