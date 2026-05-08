#ifndef SC_COMPARTMENT_H
#define SC_COMPARTMENT_H

#include <string>

namespace sc
{
    // A single physical layer in the simulation stack (vehicle or skin layer).
    // Geometry indices are filled in by Geometry::create().
    struct Compartment
    {
        std::string name;
        int    height_um   = 0;     // physical thickness, in um
        double D           = 1.0;   // um^2 / min
        double K           = 1.0;   // partition coefficient relative to vehicle
        double area_um2    = 1.0;   // cross-sectional area, in um^2
        double c_init      = 0.0;   // mg / um^3
        bool   finite_dose = true;

        int geo_from = 0;  // first space-step index belonging to this compartment
        int geo_to   = 0;  // last (inclusive)

        Compartment() = default;
        Compartment(int height, double D_, double K_, double A_, std::string n)
            : name(std::move(n)), height_um(height), D(D_), K(K_), area_um2(A_)
        {
        }
    };
}

#endif  // SC_COMPARTMENT_H
