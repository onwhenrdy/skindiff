#include "system.h"
#include "algorithms.h"
#include "helper.h"
#include <iostream>

namespace sc
{
    System::System(const Parameter& parameter)
        : m_sink_logger(parameter.systemParameter().matrixBuilderMethod(),
                        parameter.vehicleParameter().appArea() * 1.E08)
        , m_parameter(parameter)
        , m_sim_time(1)
        , m_replace_after(0)
        , m_remove_at(0)
        , m_scale(1.0)

    {
        const auto& v_params     = parameter.vehicleParameter();
        const auto& sys_params   = parameter.systemParameter();
        const auto& sink_params  = parameter.sinkParameter();
        const auto& pk_params    = parameter.pkParameter();
        const auto& layer_params = parameter.layer_parameter();
        const auto& log_params   = parameter.logParameter();

        // in minutes
        m_replace_after = v_params.replaceAfter();
        // in minutes
        m_remove_at = v_params.removeAt();
        // set sim time in minutes
        this->setSimTime(sys_params.simulationTime());

        // scaling (base is mg)
        if (log_params.scaling() == LogParameter::Scaling::UG)
        {
            m_scale = 1.0E3;
        }
        else if (log_params.scaling() == LogParameter::Scaling::NG)
        {
            m_scale = 1.0E6;
        }

        m_matrix_builder.setMaxModule(sys_params.maxModule());

        // add layers
        // in um^2 from cm^2
        const auto app_area = v_params.appArea() * 1.0E8;
        Compartment donor(v_params.height(), v_params.D(), 1.0, app_area, v_params.name());
        donor.setCInit(v_params.cInit() * 1E-12);  // im mg / um^3
        donor.setFiniteDose(v_params.finiteDose());
        this->addCompartment(donor);
        for (const auto& layer : layer_params)
        {
            Compartment comp(layer.height(), layer.D(), layer.K(), app_area * layer.crossSection(),
                             layer.name());
            comp.setCInit(layer.cInit() * 1E-12);
            this->addCompartment(comp);
        }

        Sink sink;
        sink.setType(pk_params.enabled() ? Sink::Type::PK_Compartment : Sink::Type::Perfect_Sink);
        sink.setA(app_area * (layer_params.empty() ? 1.0 : layer_params.back().crossSection()));
        sink.setVd(sink_params.Vd());
        sink.setT_half(pk_params.thalf() * 60.0);  // in minutes
        sink.setName(sink_params.name());
        sink.setCinit(sink_params.cInit() * 1E-12);
        this->setSink(sink);

        // generate Geometry
        m_geometry.setEta(sys_params.eta());
        this->createGeometry(sys_params.discMethod(), sys_params.resolution());

        // create matrices
        m_matrix_builder.setMethod(parameter.systemParameter().matrixBuilderMethod());
        m_matrix_builder.buildMatrix(compartments(), m_geometry, &m_sink);

        // init concentration vector
        this->createInitConcentrations(m_geometry, m_compartments, &m_sink);

        const auto& mass_postfix = parameter.logParameter().massFilePostfix();
        const auto& cdp_postfix  = parameter.logParameter().CDPFilePostfix();
        const auto& f_tag        = parameter.logParameter().tag();

        const auto w_dir = log_params.workingDir();
        // setup compartment logger
        m_sink_logger.setName(m_sink.name() + " Logger");
        m_sink_logger.setFilename(w_dir + f_tag + "_" + m_sink.name() + "_" + mass_postfix +
                                  ".dat");
        m_sink_logger.registerSink(&m_sink);
        m_sink_logger.setAutoLogEnabled(sink_params.log());
        m_sink_logger.setEnabled(m_sink_logger.autoLogEnabled());
        m_sink_logger.setColumn2Name("conc");
        m_sink_logger.setZip(log_params.gzipMass());
        m_sink_logger.setLogInterval(log_params.massLogInterval());

        int i = 0;
        for (auto& comp : m_compartments)
        {
            CompartmentLog2D logger(m_parameter.systemParameter().matrixBuilderMethod(),
                                    m_parameter.vehicleParameter().appArea() * 1.E8);
            logger.setName(comp.name() + " logger");
            logger.setFilename(w_dir + f_tag + "_" + comp.name() + "_" + mass_postfix + ".dat");
            logger.registerCompartment(&comp);

            const auto enabled = (i == 0) ? v_params.log() : layer_params[i - 1].log();
            logger.setAutoLogEnabled(enabled);
            logger.setEnabled(logger.autoLogEnabled());
            logger.setZip(log_params.gzipMass());
            logger.setLogInterval(log_params.massLogInterval());
            m_compartment_logger.push_back(logger);
            i++;
        }

        // setup cdp logger
        i              = 0;
        const auto& ss = m_geometry.spaceSteps();
        for (auto& comp : m_compartments)
        {
            CompartmentLog3D logger;
            logger.setName(comp.name() + " CDP logger");
            logger.setFilename(w_dir + f_tag + "_" + comp.name() + "_" + cdp_postfix + ".dat");
            logger.registerCompartment(&comp);

            const auto enabled = (i == 0) ? v_params.logCDP() : layer_params[i - 1].logCDP();
            logger.setAutoLogEnabled(enabled);
            logger.setEnabled(logger.autoLogEnabled());
            logger.setZip(log_params.gzipCDP());
            logger.setLogInterval(log_params.CDPLogInterval());
            logger.setConcentrationPosition(m_matrix_builder.method());
            logger.setStepSizes(std::vector<double>(ss.begin() + comp.geometryFromIdx(),
                                                    ss.begin() + comp.geometryToIdx() + 1));
            m_cdp_logger.push_back(logger);
            i++;
        }
    }

    const std::vector<Compartment>& System::compartments() const
    {
        return m_compartments;
    }

    void System::addCompartment(const Compartment& compartment)
    {
        m_compartments.push_back(compartment);
    }

    const Sink& System::sink() const
    {
        return m_sink;
    }

    void System::setSink(const Sink& sink)
    {
        m_sink = sink;
    }

    void System::createGeometry(Geometry::DiscMethod method, int n_ss_per_um)
    {
        m_geometry.create(method, m_compartments, n_ss_per_um, &m_sink);
    }

    void System::createInitConcentrations(const Geometry& geometry,
                                          const std::vector<Compartment>& compartments, Sink* sink)
    {
        m_concentrations.clear();
        const auto c_size = geometry.spaceSteps().size();
        m_concentrations.resize(c_size, 0.0);
        for (const auto& comp : compartments)
        {
            const auto start_idx = comp.geometryFromIdx();
            const auto stop_idx  = comp.geometryToIdx();
            const auto comp_c    = comp.cInit();

            for (int i = start_idx; i <= stop_idx; ++i)
            {
                m_concentrations[i] = comp_c;
            }
        }

        if (sink)
        {
            // in um
            const auto ss = geometry.spaceSteps()[sink->geometryFromIdx()];
            const auto A  = sink->A();            // in um^3
            const auto Vd = sink->Vd() * 1.0E12;  // in um^3
            m_concentrations[sink->geometryFromIdx()] = sink->cInit() * Vd / (ss * A);
        }
    }

    void System::resetCompartmentConcentration(const Compartment& comp,
                                               std::vector<double>& concentrations)
    {
        const auto start_idx = comp.geometryFromIdx();
        const auto stop_idx  = comp.geometryToIdx();
        const auto comp_c    = comp.cInit();
        for (int i = start_idx; i <= stop_idx; ++i)
        {
            concentrations[i] = comp_c;
        }
    }

    void System::removeTopCompartment()
    {
        const auto top_comp = m_compartments[0];

        // remove compartment
        m_compartments.erase(m_compartments.begin());

        // rewire compartment logger
        m_compartment_logger[0].registerCompartment(nullptr);
        m_cdp_logger[0].registerCompartment(nullptr);
        for (int i = 1; i < m_compartment_logger.size(); ++i)
        {
            m_compartment_logger[i].registerCompartment(&m_compartments[i - 1]);
            m_cdp_logger[i].registerCompartment(&m_compartments[i - 1]);
        }

        // adjust geometry
        m_geometry.remove(top_comp.geometryFromIdx(), top_comp.geometryToIdx() + 1);

        // adjust concentrations
        m_concentrations.erase(m_concentrations.begin(),
                               m_concentrations.begin() + top_comp.geometryToIdx() + 1);

        // adjust compartment idxs
        const auto top_comp_size = top_comp.geometryToIdx() + 1;
        for (auto& comp : m_compartments)
        {
            comp.setGeometryIdx(comp.geometryFromIdx() - top_comp_size,
                                comp.geometryToIdx() - top_comp_size);
        }
        m_sink.setGeometryIdx(m_sink.geometryFromIdx() - top_comp_size,
                              m_sink.geometryToIdx() - top_comp_size);

        // build new matrix
        m_matrix_builder.buildMatrix(compartments(), m_geometry, &m_sink);
    }

    void System::log(double time)
    {
        m_sink_logger.log(time, m_geometry, m_concentrations, m_scale);
        for (auto& logger : m_compartment_logger)
        {
            logger.log(time, m_geometry, m_concentrations, m_scale);
        }

        for (auto& logger : m_cdp_logger)
        {
            // we log in scale / ml
            logger.log(time, m_concentrations, m_scale * 1.0E12);
        }
    }

    void System::initLogger()
    {
        const auto s_time = this->simTime();
        for (auto& logger : m_compartment_logger)
        {
            logger.setTimeHint(s_time);
        }
        m_sink_logger.setTimeHint(s_time);

        for (auto& logger : m_cdp_logger)
        {
            logger.setTimeHint(s_time);
        }
    }

    bool System::initRun()
    {
        return true;
    }

    bool System::tearDownRun()
    {
        return true;
    }

    void System::progressCallback(int current_iteration)
    {
        (void)current_iteration;
        // do nothing
    }

    bool System::testForStop(int current_iteration)
    {
        (void)current_iteration;
        return false;
    }

    const std::vector<CompartmentLog2D>& System::compartmentLogger() const
    {
        return m_compartment_logger;
    }

    bool System::writeLogsToFiles() const
    {
        bool ok = true;

        // compartment logs
        if (m_sink_logger.enabled())
        {
            ok = m_sink_logger.writeToFile();
        }

        if (ok)
        {
            for (const auto& logger : m_compartment_logger)
            {
                if (logger.enabled())
                {
                    ok = logger.writeToFile();
                    if (!ok)
                    {
                        break;
                    }
                }
            }
        }

        // cdp logs
        if (ok)
        {
            for (const auto& logger : m_cdp_logger)
            {
                if (logger.enabled())
                {
                    ok = logger.writeToFile();
                    if (!ok)
                    {
                        break;
                    }
                }
            }
        }

        return ok;
    }

    const Parameter& System::parameter() const
    {
        return m_parameter;
    }

    const CompartmentLog2D& System::sinkLogger() const
    {
        return m_sink_logger;
    }

    const std::vector<CompartmentLog3D>& System::cdpLogger() const
    {
        return m_cdp_logger;
    }

    int System::simTime() const
    {
        return m_sim_time;
    }

    void System::setSimTime(int sim_time)
    {
        if (sim_time >= 0)
        {
            m_sim_time = sim_time;
        }
    }

    const std::vector<double>& System::concentrations() const
    {
        return m_concentrations;
    }

    System::Result System::run()
    {
        if (!this->initRun())
        {
            return Result::Failed;
        }
        this->initLogger();

        auto n_ts       = m_matrix_builder.timesteps();
        auto rhs_matrix = m_matrix_builder.matrixRhs();
        auto lhs_matrix = m_matrix_builder.matrixLhs();

        auto vehicle_removed    = false;
        const auto must_replace = m_replace_after != 0;
        const auto must_remove  = m_remove_at != 0;

        // main run loop
        this->log(0);

        for (int t = 1; t <= m_sim_time; ++t)
        {
            if (testForStop(t))
            {
                return Result::Stopped;
            }
            progressCallback(t);

            for (int ts = 1; ts <= n_ts; ++ts)
            {
                // rhs vector
                rhs_matrix.inlineMultiply(m_concentrations);
                algorithm::thomasReUseIP(lhs_matrix, m_concentrations);
            }

            if (must_replace && !vehicle_removed && t > 1 && t % m_replace_after == 0)
            {
                // replace the first compartment
                this->resetCompartmentConcentration(m_compartments[0], m_concentrations);
            }

            if (must_remove && t == m_remove_at)
            {
                vehicle_removed = true;
                this->removeTopCompartment();
                rhs_matrix = m_matrix_builder.matrixRhs();
                lhs_matrix = m_matrix_builder.matrixLhs();
                n_ts       = m_matrix_builder.timesteps();
            }

            // log
            this->log(t);
        }

        if (!this->tearDownRun())
        {
            return Result::Failed;
        }

        return Result::Executed;
    }

    const Geometry& System::geometry() const
    {
        return m_geometry;
    }
}
