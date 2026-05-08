#ifndef SC_SINK_H
#define SC_SINK_H

#include <string>

namespace sc
{
    class Sink
    {
      public:
        enum class Type
        {
            Perfect_Sink,
            PK_Compartment
        };

        Sink();
        Sink(Type type, double A, double Vd, double t_half, const std::string& name);

        Type type() const;
        void setType(const Type& type);

        // in ml
        double Vd() const;
        void setVd(double Vd);

        // in min
        double t_half() const;
        void setT_half(double t_half);

        // in 1/min
        double kEl() const;
        void setkEl(double val);

        std::string name() const;
        void setName(const std::string& name);

        // in um^2
        double A() const;
        void setA(double A);

        int geometryFromIdx() const;
        int geometryToIdx() const;
        void setGeometryIdx(int from , int to);

        // in mg/um^3
        double cInit() const;
        void setCinit(double c_init);

    private:
        std::string m_name;
        Type m_type;
        double m_A;       // in um^2
        double m_Vd;      // in ml
        double m_t_half;  // in min
        int m_geo_from;
        int m_geo_to;
        double m_c_init;  // in mg/ml
    };
}
#endif  // SINK_H
