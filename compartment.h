#ifndef SC_COMPARTMENT_H
#define SC_COMPARTMENT_H

#include <string>

namespace sc
{
    class Compartment
    {
      public:
        Compartment();
        Compartment(int size, double D, double K, double A, const std::string& name);

        std::string name() const;
        void setName(const std::string& name);

        // in um
        int size() const;
        void setSize(int size);

        // in um^2/min
        double D() const;
        void setD(double D);

        // in um^2
        double A() const;
        void setA(double A);

        int geometryFromIdx() const;
        int geometryToIdx() const;
        void setGeometryIdx(int from , int to);

        // in mg/um^3
        double cInit() const;
        void setCInit(double value);

        double K() const;
        void setK(double K);

        bool finiteDose() const;
        void setFiniteDose(bool finite_dose);

    private:
        std::string m_name;
        int m_size;  // in um
        double m_D;  // in um^2/min
        double m_K;
        double m_A;  // in um^2
        int m_geo_from;
        int m_geo_to;
        double m_c_init; // in mg/um^3
        bool m_finite_dose;
    };
}
#endif  // COMPARTMENT_H
