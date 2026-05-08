#ifndef SC_SINK_H
#define SC_SINK_H

#include <cmath>
#include <string>

namespace sc
{
    // The bottom boundary of the compartment stack.
    // Geometry indices are filled in by Geometry::create().
    struct Sink
    {
        enum class Type
        {
            Perfect_Sink,
            PK_Compartment
        };

        std::string name;
        Type   type     = Type::Perfect_Sink;
        double area_um2 = 1.0;   // um^2
        double Vd       = 1.0;   // ml
        double t_half   = 1.0;   // min
        double c_init   = 0.0;   // mg / um^3

        int geo_from = 0;
        int geo_to   = 0;

        // Elimination rate constant (1 / min)
        [[nodiscard]] double kEl() const noexcept
        {
            return std::log(2.0) / t_half;
        }
    };
}

#endif  // SC_SINK_H
