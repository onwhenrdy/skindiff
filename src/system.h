#ifndef SC_SYSTEM_H
#define SC_SYSTEM_H

#include "compartment.h"
#include "geometry.h"
#include "matrixbuilder.h"
#include "parameter.h"
#include "sink.h"
#include "compartmentlog2d.h"
#include "compartmentlog3d.h"
#include <vector>

namespace sc
{
    class System
    {
      public:
        enum class Status
        {
            Idle,
            Runs
        };

        enum class Result
        {
            Executed,
            Failed,
            Stopped
        };


        System(const Parameter& parameter);
        virtual ~System() = default;

        const std::vector<Compartment>& compartments() const;
        const Sink& sink() const;
        const Geometry& geometry() const;
        const std::vector<double>& concentrations() const;

        Result run();

        int simTime() const;
        void setSimTime(int sim_time);

        const CompartmentLog2D& sinkLogger() const;
        const std::vector<CompartmentLog2D>& compartmentLogger() const;
        const std::vector<CompartmentLog3D>& cdpLogger() const;
        bool writeLogsToFiles() const;
        const Parameter& parameter() const;

      private:
        void addCompartment(const Compartment& compartment);
        void setSink(const Sink& sink);

        void createGeometry(Geometry::DiscMethod method, int n_ss_per_um);
        void createInitConcentrations(const Geometry& geometry,
                                      const std::vector<Compartment>& compartments,
                                      Sink* sink = nullptr);

        void resetCompartmentConcentration(const Compartment& comp,
                                           std::vector<double>& concentrations);
        void removeTopCompartment();

        void log(double time);
        void initLogger();

      protected:
        virtual bool initRun();
        virtual bool tearDownRun();
        virtual void progressCallback(int current_iteration);
        virtual bool testForStop(int current_iteration);

      private:
        // data
        std::vector<Compartment> m_compartments;
        std::vector<double> m_concentrations;
        Sink m_sink;
        Geometry m_geometry;
        MatrixBuilder m_matrix_builder;

        // logger
        std::vector<CompartmentLog2D> m_compartment_logger;
        CompartmentLog2D m_sink_logger;
        std::vector<CompartmentLog3D> m_cdp_logger;

        // properties
        Parameter m_parameter;
        int m_sim_time;       // in min
        int m_replace_after;  // in min
        int m_remove_at;      // in min
        double m_scale;       // sink scaling factor
    };
}
#endif  // SYSTEM_H
