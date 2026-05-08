#include "compartment.h"
#include <cassert>

namespace sc
{
    Compartment::Compartment()
        : m_size(0)
        , m_D(1.0)
        , m_K(1.0)
        , m_A(1.0)
        , m_geo_from(0)
        , m_geo_to(0)
        , m_c_init(0.0)
        , m_finite_dose(true)
    {
    }

    Compartment::Compartment(int size, double D, double K, double A, const std::string& name)
        : m_name(name)
        , m_size(size)
        , m_D(D)
        , m_K(K)
        , m_A(A)
        , m_geo_from(0)
        , m_geo_to(0)
        , m_c_init(0.0)
        , m_finite_dose(true)
    {
    }

    std::string Compartment::name() const
    {
        return m_name;
    }

    void Compartment::setName(const std::string& name)
    {
        m_name = name;
    }

    int Compartment::size() const
    {
        return m_size;
    }

    void Compartment::setSize(int size)
    {
        if (size >= 1)
        {
            m_size = size;
        }
    }

    double Compartment::D() const
    {
        return m_D;
    }

    void Compartment::setD(double D)
    {
        if (D >= 0)
        {
            m_D = D;
        }
    }

    double Compartment::A() const
    {
        return m_A;
    }

    void Compartment::setA(double A)
    {
        if (A >= 0)
        {
            m_A = A;
        }
    }

    int Compartment::geometryFromIdx() const
    {
        return m_geo_from;
    }

    int Compartment::geometryToIdx() const
    {
        return m_geo_to;
    }

    void Compartment::setGeometryIdx(int from, int to)
    {
        assert(from <= to);
        m_geo_from = from;
        m_geo_to   = to;
    }

    double Compartment::cInit() const
    {
        return m_c_init;
    }

    void Compartment::setCInit(double value)
    {
        assert(value >= 0.0);
        m_c_init = value;
    }

    double Compartment::K() const
    {
        return m_K;
    }

    void Compartment::setK(double K)
    {
        if (K > 0)
        {
            m_K = K;
        }
    }

    bool Compartment::finiteDose() const
    {
        return m_finite_dose;
    }

    void Compartment::setFiniteDose(bool finite_dose)
    {
        m_finite_dose = finite_dose;
    }
}
