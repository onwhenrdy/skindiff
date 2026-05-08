#include "sink.h"
#include <cassert>
#include <cmath>

namespace sc
{
    Sink::Sink()
        : m_type(Sink::Type::Perfect_Sink)
        , m_A(1.0)
        , m_Vd(1.0)
        , m_t_half(1.0)
        , m_geo_from(0)
        , m_geo_to(0)
        , m_c_init(0.0)
    {
    }

    Sink::Sink(Sink::Type type, double A, double Vd, double t_half, const std::string& name)
        : m_name(name)
        , m_type(type)
        , m_A(A)
        , m_Vd(Vd)
        , m_t_half(t_half)
        , m_geo_from(0)
        , m_geo_to(0)
        , m_c_init(0.0)
    {
    }

    Sink::Type Sink::type() const
    {
        return m_type;
    }

    void Sink::setType(const Sink::Type& type)
    {
        m_type = type;
    }

    double Sink::Vd() const
    {
        return m_Vd;
    }

    void Sink::setVd(double Vd)
    {
        if (Vd > 0)
        {
            m_Vd = Vd;
        }
    }

    double Sink::t_half() const
    {
        return m_t_half;
    }

    void Sink::setT_half(double t_half)
    {
        if (t_half > 0)
        {
            m_t_half = t_half;
        }
    }

    double Sink::kEl() const
    {
        return std::log(2) / m_t_half;
    }

    void Sink::setkEl(double val)
    {
        if (val > 0.0)
        {
            m_t_half = std::log(2) / val;
        }
    }

    std::string Sink::name() const
    {
        return m_name;
    }

    void Sink::setName(const std::string& name)
    {
        m_name = name;
    }

    double Sink::A() const
    {
        return m_A;
    }

    void Sink::setA(double A)
    {
        if (A > 0)
        {
            m_A = A;
        }
    }

    int Sink::geometryFromIdx() const
    {
        return m_geo_from;
    }

    int Sink::geometryToIdx() const
    {
        return m_geo_to;
    }

    void Sink::setGeometryIdx(int from, int to)
    {
        assert(from <= to);
        m_geo_from = from;
        m_geo_to   = to;
    }

    double Sink::cInit() const
    {
        return m_c_init;
    }

    void Sink::setCinit(double c_init)
    {
        m_c_init = c_init;
    }
}
