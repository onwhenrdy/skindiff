#ifndef SC_SINK_H
#define SC_SINK_H

#include <string>

namespace sc
{
    // The bottom boundary of the compartment stack.
    //
    // The sink behaves as a perfect Dirichlet boundary for the diffusion
    // problem -- the membrane sees u_sink = 0 regardless of `Vd`. The sink
    // cell itself is a virtual mass accumulator; its stored value is
    // Vd-scaled so that integrating it recovers the cumulative mass that
    // crossed the membrane.
    //
    // Geometry indices are filled in by Geometry::create().
    struct Sink
    {
        std::string name;
        double area_um2 = 1.0;   // um^2
        double Vd       = 1.0;   // ml
        double c_init   = 0.0;   // mg / um^3

        int geo_from = 0;
        int geo_to   = 0;
    };
}

#endif  // SC_SINK_H
