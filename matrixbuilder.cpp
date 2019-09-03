#include "matrixbuilder.h"
#include <algorithm>
#include <cassert>
#include <cmath>

namespace sc
{
    MatrixBuilder::MatrixBuilder(Method method)
        : m_method(method), m_max_module(50.0), m_timesteps(1)
    {
    }

    MatrixBuilder::Method MatrixBuilder::method() const
    {
        return m_method;
    }

    void MatrixBuilder::setMethod(const MatrixBuilder::Method& method)
    {
        m_method = method;
    }

    double MatrixBuilder::maxModule() const
    {
        return m_max_module;
    }

    void MatrixBuilder::setMaxModule(double max_module)
    {
        assert(max_module > 0.0);
        m_max_module = max_module;
    }

    const TDMatrix& MatrixBuilder::matrixRhs() const
    {
        return m_matrix_rhs;
    }

    const TDMatrix& MatrixBuilder::matrixLhs() const
    {
        return m_matrix_lhs;
    }

    bool MatrixBuilder::buildMatrix(const std::vector<Compartment>& compartments,
                                    const Geometry& geometry, Sink* sink)
    {
        switch (m_method)
        {
            case MatrixBuilder::Method::DSkin_1_3:
                return this->build_M_DS_1_3(compartments, geometry, sink);
            case MatrixBuilder::Method::DSkin_1_4:
                return this->build_M_DS_1_4(compartments, geometry, sink);
            case MatrixBuilder::Method::DSkin_1_5:
                return this->build_M_DS_1_5(compartments, geometry, sink);
            default:
                assert(false);
                break;
        }

        return false;
    }

    int MatrixBuilder::timesteps() const
    {
        return m_timesteps;
    }

    bool MatrixBuilder::build_M_DS_1_3(const std::vector<Compartment>& compartments,
                                       const Geometry& geometry, Sink* sink)
    {
        assert(compartments.size() > 0);

        const auto sys_size = geometry.size();

        const auto D_vec = this->createParamVector(
            sys_size, compartments, [](const Compartment& comp) { return comp.D(); }, sink);
        const auto K_vec = this->createParamVector(
            sys_size, compartments, [](const Compartment& comp) { return comp.K(); }, sink);
        const auto A_vec = this->createParamVector(
            sys_size, compartments, [](const Compartment& comp) { return comp.A(); }, sink);

        m_matrix_rhs = TDMatrix(sys_size);
        // last diag idx: size - 1
        // last upper/lower idx: size - 2
        // -> sink

        const auto space_steps = geometry.spaceSteps();

        // -------------------------------------------------------------
        // BUILD MAIN MATRIX ELEMENTS
        //
        // Matrix is build as "Element_i = 2.0 * M_i"
        //
        // with M = D * dt/(dx)^2
        //
        // first elements at the TOP (boundary conditions)
        // -> reflecting boundaries (flux is 0)
        // dc/dx_(x=0) = 0
        auto l_dx            = space_steps[0];
        auto r_dx            = this->avgFromIdx(space_steps, 0, 1);
        m_matrix_rhs.diag(0) = 2.0 * D_vec[0] / (l_dx * r_dx);
        // gain 2 times from right to conserve the mass since we spend
        // mass to the bounday in the upper equation
        m_matrix_rhs.upper(0) = D_vec[0] * 4.0 / (r_dx * (l_dx + r_dx));

        // main building loop for finite dose
        for (int i = 1; i < sys_size - 1; ++i)
        {
            auto l   = this->avgFromIdx(space_steps, i, i - 1);
            auto r   = this->avgFromIdx(space_steps, i, i + 1);
            auto D_r = this->harmMeanFromIdx(D_vec, i, i + 1);
            auto D_l = this->harmMeanFromIdx(D_vec, i, i - 1);

            // check for backflux reduction
            auto k1 = 1.0;
            auto k2 = 1.0;
            auto k3 = 1.0;
            auto k4 = 1.0;
            this->backFluxCorrection(K_vec, i, k1, k2, k3, k4);

            // area correction
            auto v1 = 1.0;
            auto v2 = 1.0;
            this->areaCorrection(A_vec, i, v1, v2);

            const auto lower_val = D_l * k1 * v2 * 2.0 / (l * (l + r));
            const auto mid_val   = (D_l * k3 * v2 + D_r * k4 * v1) / (l * r);
            const auto upper_val = D_r * k2 * v1 * 2.0 / (r * (l + r));

            m_matrix_rhs.diag(i)      = mid_val;
            m_matrix_rhs.upper(i)     = upper_val;
            m_matrix_rhs.lower(i - 1) = lower_val;
        }

        // last matrix row element for the lower part (no boundaries applied yet !!)
        l_dx     = this->avgFromIdx(space_steps, sys_size - 1, sys_size - 2);
        r_dx     = space_steps[sys_size - 1];
        auto D_l = D_vec[sys_size - 1];
        m_matrix_rhs.lower(sys_size - 2) = D_l * 2.0 / (l_dx * (l_dx + r_dx));

        // -------------------------------------------------------------
        // GET MAX module, get ts and apply dt
        const auto max_m = m_matrix_rhs.absMax();
        // assign nr of ts
        m_timesteps   = static_cast<int>(std::ceil(max_m / m_max_module));
        const auto dt = 1.0 / m_timesteps;
        // include dt to all
        m_matrix_rhs.multiplyBy(dt);

        // correct the matrix equations
        for (int i = 0; i < sys_size - 1; ++i)
        {
            m_matrix_rhs.diag(i)  = 1.0 - m_matrix_rhs.diag(i) / 2.0;
            m_matrix_rhs.lower(i) = m_matrix_rhs.lower(i) / 2.0;
            m_matrix_rhs.upper(i) = m_matrix_rhs.upper(i) / 2.0;
        }
        // -------------------------------------------------------------
        // Lower boundary conditions to the sink element
        // do not get mass from the last element (only spend it)
        m_matrix_rhs.upper(sys_size - 2) = 0.0;

        // -------------------------------------------------------------
        // Handle sink element and kinetics
        if (sink)
        {
            // conserve the the concentration in the sink
            m_matrix_rhs.diag(sys_size - 1) = 1.0;
            if (sink->type() == Sink::Type::PK_Compartment)
            {
                m_matrix_rhs.diag(sys_size - 1) = (1.0 - dt * sink->kEl() / 2.0);
            }
        }

        // Infinite dose settings
        if (!compartments[0].finiteDose())
        {
            assert(false);
        }

        // BUILD LHS matrix
        m_matrix_lhs = this->fromRhs(m_matrix_rhs);
        return true;
    }

    bool MatrixBuilder::build_M_DS_1_4(const std::vector<Compartment>& compartments,
                                       const Geometry& geometry, Sink* sink)
    {
        assert(compartments.size() > 0);

        const auto sys_size = geometry.size();

        const auto D_vec = this->createParamVector(
            sys_size, compartments, [](const Compartment& comp) { return comp.D(); }, sink);
        const auto K_vec = this->createParamVector(
            sys_size, compartments, [](const Compartment& comp) { return comp.K(); }, sink);
        const auto A_vec = this->createParamVector(
            sys_size, compartments, [](const Compartment& comp) { return comp.A(); }, sink);

        m_matrix_rhs = TDMatrix(sys_size);
        // last diag idx: size - 1
        // last upper/lower idx: size - 2
        // -> sink

        const auto space_steps = geometry.spaceSteps();

        m_matrix_lhs = TDMatrix(sys_size);

        // Reflecting boundary at the start
        auto l_c = space_steps[0];
        auto l_r = space_steps[1];

        auto D_c = D_vec[0];
        auto D_r = D_vec[1];

        auto K_c = K_vec[0];
        auto K_r = K_vec[1];

        auto A_c = A_vec[0];
        auto A_r = A_vec[1];

        auto h2      = (l_c + l_r) / 2.0;
        auto upper_f = (l_c + l_r) * D_c * D_r / (l_c * D_r + K_c / K_r * l_r * D_c) /
                       (h2 * h2);  // subsequent element

        // Gain from subsequent element
        auto upper_val = upper_f * K_c / K_r * std::min(1.0, A_r / A_c);

        // Loss of central element
        auto mid_val = upper_f * std::min(1.0, A_r / A_c);

        m_matrix_rhs.diag(0)  = mid_val;
        m_matrix_rhs.upper(0) = upper_val;

        // Main building loop for finite dose
        for (int i = 1; i < sys_size - 1; ++i)
        {
            // Properties of left, central and right element
            auto l_l = space_steps[i - 1];
            auto l_c = space_steps[i];
            auto l_r = space_steps[i + 1];

            auto D_l = D_vec[i - 1];
            auto D_c = D_vec[i];
            auto D_r = D_vec[i + 1];

            auto K_l = K_vec[i - 1];
            auto K_c = K_vec[i];
            auto K_r = K_vec[i + 1];

            auto A_l = A_vec[i - 1];
            auto A_c = A_vec[i];
            auto A_r = A_vec[i + 1];

            auto h1 = (l_l + l_c) / 2.0;
            auto h2 = (l_c + l_r) / 2.0;

            // The calculation uses the weighted harmonic mean adjusted for the partition coeffient
            // times the correction for the space steps
            // The weighted harmonic mean computed below is only valid for equal space steps, thus
            // we need 4 equally long space steps at the boundary!
            auto lower_f = (l_l + l_c) * D_l * D_c / (l_l * D_c + K_l / K_c * l_c * D_l) * 2.0 *
                           h2 / (h1 * h2 * (h1 + h2));  // preceding element
            auto upper_f = (l_c + l_r) * D_c * D_r / (l_c * D_r + K_c / K_r * l_r * D_c) * 2.0 *
                           h1 / (h1 * h2 * (h1 + h2));  // subsequent element

            auto lower_val = lower_f * std::min(1.0, A_l / A_c);
            auto upper_val = upper_f * K_c / K_r * std::min(1.0, A_r / A_c);
            auto mid_val =
                lower_f * K_l / K_c * std::min(1.0, A_l / A_c) + upper_f * std::min(1.0, A_r / A_c);

            m_matrix_rhs.diag(i)      = mid_val;
            m_matrix_rhs.upper(i)     = upper_val;
            m_matrix_rhs.lower(i - 1) = lower_val;
        }

        // Last matrix row element for the lower part (no boundaries applied yet !!)
        auto l_l = space_steps[sys_size - 2];
        l_c      = space_steps[sys_size - 1];

        auto D_l = D_vec[sys_size - 2];
        D_c      = D_vec[sys_size - 1];

        auto K_l = K_vec[sys_size - 2];
        K_c      = K_vec[sys_size - 1];

        auto A_l = A_vec[sys_size - 2];
        A_c      = A_vec[sys_size - 1];

        // Define helper variables

        auto h1      = (l_l + l_c) / 2.0;
        auto lower_f = (l_l + l_c) * D_l * D_c / (l_l * D_c + K_l / K_c * l_c * D_l) /
                       (h1 * h1);  // preceding element

        // Gain from preceding element
        auto lower_val = lower_f * std::min(1.0, A_l / A_c);

        // Loss of central element
        mid_val = lower_f * K_l / K_c * std::min(1.0, A_l / A_c);

        m_matrix_rhs.diag(sys_size - 1)  = mid_val;  // reflecting boundary
        m_matrix_rhs.lower(sys_size - 2) = lower_val;

        // -------------------------------------------------------------
        // Get MAX module, compute ts and apply dt
        const auto max_m = m_matrix_rhs.absMax();

        // Compute number of time steps per minute
        m_timesteps = static_cast<int>(std::max(1.0, std::ceil(max_m / m_max_module)));

        const auto dt = 1.0 / m_timesteps;

        // Mutiply preliminary matrix by dt
        m_matrix_rhs.multiplyBy(dt);

        // Correct the matrix and populate lhs matrix
        for (int i = 0; i < sys_size - 1; ++i)
        {
            m_matrix_lhs.diag(i)  = 2.0 + m_matrix_rhs.diag(i);
            m_matrix_lhs.lower(i) = -m_matrix_rhs.lower(i);
            m_matrix_lhs.upper(i) = -m_matrix_rhs.upper(i);

            m_matrix_rhs.diag(i) = 2.0 - m_matrix_rhs.diag(i);
        }

        // -------------------------------------------------------------
        // Lower boundary conditions to the sink element
        // do not gain mass from the last element (only spend to it)
        m_matrix_rhs.upper(sys_size - 2) = 0.0;
        m_matrix_lhs.upper(sys_size - 2) = 0.0;

        // -------------------------------------------------------------
        // Handle sink element and kinetics
        if (sink)
        {
            // conserve the concentration in the sink

            m_matrix_rhs.diag(sys_size - 1) = 2.0;
            m_matrix_lhs.diag(sys_size - 1) = 2.0;

            if (sink->type() == Sink::Type::PK_Compartment)
            {
                m_matrix_rhs.diag(sys_size - 1) = 2.0 - dt * sink->kEl();
                m_matrix_lhs.diag(sys_size - 1) = 2.0 + dt * sink->kEl();
            }
        }

        // Infinite dose settings
        if (!compartments[0].finiteDose())
        {
            m_matrix_rhs.diag(0) = 2.0;
            m_matrix_lhs.diag(0) = 2.0;

            m_matrix_rhs.upper(0) = 0.0;
            m_matrix_lhs.upper(0) = 0.0;
        }

        return true;
    }

    bool MatrixBuilder::build_M_DS_1_5(const std::vector<Compartment>& compartments,
                                       const Geometry& geometry, Sink* sink)
    {
        assert(compartments.size() > 0);

        const auto sys_size = geometry.size();

        const auto D_vec = this->createParamVector(
            sys_size, compartments, [](const Compartment& comp) { return comp.D(); }, sink);
        const auto K_vec = this->createParamVector(
            sys_size, compartments, [](const Compartment& comp) { return comp.K(); }, sink);
        const auto A_vec = this->createParamVector(
            sys_size, compartments, [](const Compartment& comp) { return comp.A(); }, sink);

        m_matrix_rhs = TDMatrix(sys_size);
        // last diag idx: size - 1
        // last upper/lower idx: size - 2
        // -> sink

        const auto space_steps = geometry.spaceSteps();

        // Code corrsponding to Crank

        m_matrix_lhs = TDMatrix(sys_size);

        // Reflecting boundary at the start
        auto l_c = space_steps[0];
        auto l_r = space_steps[1];

        auto D_c = D_vec[0];
        auto D_r = D_vec[1];

        auto K_c = K_vec[0] * A_vec[0];
        auto K_r = K_vec[1] * A_vec[1];

        auto h2      = (l_c + l_r) / 2.0;
        auto upper_f = (l_c + l_r) * D_c * D_r / (l_c * D_r + K_c / K_r * l_r * D_c) /
                       (h2 * h2);  // subsequent element

        // Gain from subsequent element
        auto upper_val = upper_f * K_c / K_r;

        // Loss of central element
        auto mid_val = upper_f;

        m_matrix_rhs.diag(0)  = mid_val;
        m_matrix_rhs.upper(0) = upper_val;

        // Main building loop for finite dose
        for (int i = 1; i < sys_size - 1; ++i)
        {
            // Properties of left, central and right element
            auto l_l = space_steps[i - 1];
            auto l_c = space_steps[i];
            auto l_r = space_steps[i + 1];

            auto D_l = D_vec[i - 1];
            auto D_c = D_vec[i];
            auto D_r = D_vec[i + 1];

            auto K_l = K_vec[i - 1] * A_vec[i - 1];
            auto K_c = K_vec[i] * A_vec[i];
            auto K_r = K_vec[i + 1] * A_vec[i + 1];
            auto h1  = (l_l + l_c) / 2.0;
            auto h2  = (l_c + l_r) / 2.0;

            // The calculation uses the weighted harmonic mean adjusted for the partition coeffient
            // times the correction for the space steps
            // The weighted harmonic mean computed below is only valid for equal space steps, thus
            // we need 4 equally long space steps at the boundary!
            auto lower_f = (l_l + l_c) * D_l * D_c / (l_l * D_c + K_l / K_c * l_c * D_l) * 2.0 *
                           h2 / (h1 * h2 * (h1 + h2));  // preceding element
            auto upper_f = (l_c + l_r) * D_c * D_r / (l_c * D_r + K_c / K_r * l_r * D_c) * 2.0 *
                           h1 / (h1 * h2 * (h1 + h2));  // subsequent element

            auto lower_val            = lower_f;
            auto upper_val            = upper_f * K_c / K_r;
            auto mid_val              = lower_f * K_l / K_c + upper_f;
            m_matrix_rhs.diag(i)      = mid_val;
            m_matrix_rhs.upper(i)     = upper_val;
            m_matrix_rhs.lower(i - 1) = lower_val;
        }

        // Last matrix row element for the lower part (no boundaries applied yet !!)
        auto l_l = space_steps[sys_size - 2];
        l_c      = space_steps[sys_size - 1];

        auto D_l = D_vec[sys_size - 2];
        D_c      = D_vec[sys_size - 1];

        auto K_l = K_vec[sys_size - 2] * A_vec[sys_size - 2];
        K_c      = K_vec[sys_size - 1] * A_vec[sys_size - 1];

        // Define helper variables
        auto h1      = (l_l + l_c) / 2.0;
        auto lower_f = (l_l + l_c) * D_l * D_c / (l_l * D_c + K_l / K_c * l_c * D_l) /
                       (h1 * h1);  // preceding element

        // Gain from preceding element
        auto lower_val = lower_f;

        // Loss of central element
        mid_val = lower_f * K_l / K_c;

        m_matrix_rhs.diag(sys_size - 1)  = mid_val;  // reflecting boundary
        m_matrix_rhs.lower(sys_size - 2) = lower_val;

        // -------------------------------------------------------------
        // Get MAX module, compute ts and apply dt
        const auto max_m = m_matrix_rhs.absMax();

        // Compute number of time steps per minute
        m_timesteps = static_cast<int>(std::max(1.0, std::ceil(max_m / m_max_module)));

        const auto dt = 1.0 / m_timesteps;

        // Mutiply preliminary matrix by dt
        m_matrix_rhs.multiplyBy(dt);

        // Correct the matrix and populate lhs matrix
        for (int i = 0; i < sys_size - 1; ++i)
        {
            m_matrix_lhs.diag(i)  = 2.0 + m_matrix_rhs.diag(i);
            m_matrix_lhs.lower(i) = -m_matrix_rhs.lower(i);
            m_matrix_lhs.upper(i) = -m_matrix_rhs.upper(i);

            m_matrix_rhs.diag(i) = 2.0 - m_matrix_rhs.diag(i);
        }

        // -------------------------------------------------------------
        // Lower boundary conditions to the sink element
        // do not gain mass from the last element (only spend to it)
        m_matrix_rhs.upper(sys_size - 2) = 0.0;
        m_matrix_lhs.upper(sys_size - 2) = 0.0;

        // -------------------------------------------------------------
        // Handle sink element and kinetics
        if (sink)
        {
            // conserve the concentration in the sink

            m_matrix_rhs.diag(sys_size - 1) = 2.0;
            m_matrix_lhs.diag(sys_size - 1) = 2.0;

            if (sink->type() == Sink::Type::PK_Compartment)
            {
                m_matrix_rhs.diag(sys_size - 1) = 2.0 - dt * sink->kEl();
                m_matrix_lhs.diag(sys_size - 1) = 2.0 + dt * sink->kEl();
            }
        }

        // Infinite dose settings
        if (!compartments[0].finiteDose())
        {
            m_matrix_rhs.diag(0) = 2.0;
            m_matrix_lhs.diag(0) = 2.0;

            m_matrix_rhs.upper(0) = 0.0;
            m_matrix_lhs.upper(0) = 0.0;
        }

        return true;
    }

    std::vector<double>
    MatrixBuilder::createParamVector(int size, const std::vector<Compartment>& compartments,
                                     std::function<double(const Compartment&)> fun, Sink* sink)
    {
        std::vector<double> result(size, 0.0);

        for (const auto& comp : compartments)
        {
            const auto start_idx = comp.geometryFromIdx();
            const auto end_idx   = comp.geometryToIdx();
            const auto val       = fun(comp);

            for (int i = start_idx; i <= end_idx; ++i)
            {
                result[i] = val;
            }
        }

        if (sink)
        {
            const auto idx = sink->geometryFromIdx();
            result[idx]    = result[idx - 1];
        }

        return result;
    }

    TDMatrix MatrixBuilder::fromRhs(const TDMatrix& rhs)
    {
        auto res        = rhs;
        const auto size = res.size();
        for (int i = 0; i < size - 1; ++i)
        {
            res.diag(i)  = 2.0 - res.diag(i);
            res.upper(i) = -res.upper(i);
            res.lower(i) = -res.lower(i);
        }
        res.diag(size - 1) = 2.0 - res.diag(size - 1);

        return res;
    }

    std::string toString(MatrixBuilder::Method method)
    {
        switch (method)
        {
            case MatrixBuilder::Method::DSkin_1_3:
                return "DSkin_1_3";
            case MatrixBuilder::Method::DSkin_1_4:
                return "DSkin_1_4";
            case MatrixBuilder::Method::DSkin_1_5:
                return "DSkin_1_5";
            default:
                return "unknown";
                break;
        }
    }

    MatrixBuilder::Method mbMethodFromString(const std::string& str, bool* ok)
    {
        if (ok)
        {
            *ok = true;
        }

        if (str == "DSkin_1_3")
        {
            return MatrixBuilder::Method::DSkin_1_3;
        }

        if (str == "DSkin_1_4")
        {
            return MatrixBuilder::Method::DSkin_1_4;
        }

        if (str == "DSkin_1_5")
        {
            return MatrixBuilder::Method::DSkin_1_5;
        }

        if (ok)
        {
            *ok = false;
        }

        return MatrixBuilder::Method::DSkin_1_3;
    }
}
