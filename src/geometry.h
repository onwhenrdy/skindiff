#ifndef SC_GEOMETRY_H
#define SC_GEOMETRY_H

#include "compartment.h"
#include "sink.h"
#include <vector>

namespace sc
{
    class Geometry
    {
      public:
        enum class DiscMethod
        {
            Undefined,
            EQUI_DIST,  //!< Equidistant mesh
            B_AND_K     //!< Grid refinement by Babucke & Kloker, 2009
        };

        Geometry();
        bool create(DiscMethod method, std::vector<Compartment>&, int ss_per_um,
                    Sink* sink = nullptr);
        const std::vector<double>& spaceSteps() const;
        int size() const;
        std::string spaceStepsR(const std::string& var_name = "a") const;
        double minSpaceStep() const;
        double maxSpaceStep() const;
        Geometry::DiscMethod discMethod() const;
        bool valid() const;
        void remove(int from_idx, int to_idx);

        double eta() const;
        void setEta(double eta);
        double calculatedEta() const;

    private:
        void findOptTransition(int& n, double& x, int& a, double& delta_x, double err);
        double findOptimalX(double start_x, int n, double a, double err);

        //! Computes the power series sum = sum_i=1^(n-1) (x^n) + x^n
        static double powerSeriesDoubleLastElement(int n, double x);

      private:
        std::vector<double> m_space_steps;
        double m_min_space_step;
        double m_max_space_step;
        DiscMethod m_disc_method;

        // properties
        bool m_valid;
        double m_eta;
        double m_calcuated_eta;
    };

    std::string toString(Geometry::DiscMethod method);
    Geometry::DiscMethod discMethodFromString(const std::string& str, bool* ok = nullptr);
}

#endif  // GEMOMETRY_H
