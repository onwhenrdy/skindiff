#include "parameter.h"
#include "system.h"

#include <Rcpp.h>

#include <string>
#include <utility>
#include <vector>

using namespace sc;

namespace
{
    // ---------- R list -> Parameters helpers ----------

    template <typename T>
    T pick(const Rcpp::List& list, const char* key, T fallback)
    {
        if (!list.containsElementNamed(key)) return fallback;
        SEXP s = list[key];
        if (Rf_isNull(s)) return fallback;
        return Rcpp::as<T>(s);
    }

    Scaling parseScaling(const std::string& s)
    {
        const auto v = scalingFromString(s);
        if (!v) Rcpp::stop("Unknown scaling '" + s + "' (expected 'mg', 'ug', or 'ng')");
        return *v;
    }

    SystemParams readSys(const Rcpp::List& sys)
    {
        SystemParams out;
        out.resolution      = pick<int>(sys,    "resolution",      1);
        out.max_module      = pick<double>(sys, "max_module",      50.0);
        out.simulation_time = pick<int>(sys,    "simulation_time", 600);
        return out;
    }

    LogParams readLog(const Rcpp::List& log)
    {
        LogParams out;
        out.scaling           = parseScaling(pick<std::string>(log, "scaling", "mg"));
        out.mass_log_interval = pick<int>(log, "mass_log_interval", 1);
        out.cdp_log_interval  = pick<int>(log, "cdp_log_interval",  1);
        return out;
    }

    SinkParams readSink(const Rcpp::List& s)
    {
        SinkParams out;
        out.name     = pick<std::string>(s, "name",     "Sink");
        out.c_init   = pick<double>(s,      "c_init",   0.0);
        out.Vd       = pick<double>(s,      "Vd",       1.0);
        out.log_mass = pick<bool>(s,        "log_mass", true);
        return out;
    }

    VehicleParams readVehicle(const Rcpp::List& v)
    {
        VehicleParams out;
        out.name          = pick<std::string>(v, "name",          "Vehicle");
        out.c_init        = pick<double>(v,      "c_init",        1.0);
        out.app_area      = pick<double>(v,      "app_area",      1.0);
        out.D             = pick<double>(v,      "D",             1.0);
        out.height        = pick<int>(v,         "height",        10);
        out.replace_after = pick<int>(v,         "replace_after", 0);
        out.remove_at     = pick<int>(v,         "remove_at",     0);
        out.finite_dose   = pick<bool>(v,        "finite_dose",   true);
        out.log_mass      = pick<bool>(v,        "log_mass",      true);
        out.log_cdp       = pick<bool>(v,        "log_cdp",       false);
        return out;
    }

    std::vector<LayerParams> readLayers(const Rcpp::List& layers)
    {
        std::vector<LayerParams> out;
        out.reserve(static_cast<std::size_t>(layers.size()));
        for (R_xlen_t i = 0; i < layers.size(); ++i)
        {
            const Rcpp::List l(layers[i]);
            LayerParams p;
            p.name          = pick<std::string>(l, "name",          "Layer");
            p.c_init        = pick<double>(l,      "c_init",        0.0);
            p.D             = pick<double>(l,      "D",             1.0);
            p.K             = pick<double>(l,      "K",             1.0);
            p.cross_section = pick<double>(l,      "cross_section", 1.0);
            p.height        = pick<int>(l,         "height",        10);
            p.log_mass      = pick<bool>(l,        "log_mass",      true);
            p.log_cdp       = pick<bool>(l,        "log_cdp",       false);
            out.push_back(std::move(p));
        }
        return out;
    }

    Parameters parametersFromR(const Rcpp::List& p)
    {
        Parameters out;
        if (p.containsElementNamed("sys"))     out.sys     = readSys(p["sys"]);
        if (p.containsElementNamed("log"))     out.log     = readLog(p["log"]);
        if (p.containsElementNamed("sink"))    out.sink    = readSink(p["sink"]);
        if (p.containsElementNamed("vehicle")) out.vehicle = readVehicle(p["vehicle"]);
        if (p.containsElementNamed("layers"))  out.layers  = readLayers(p["layers"]);
        return out;
    }

    // ---------- System -> R ----------

    // Wraps the simulation core to forward progress / interrupt to R.
    class SystemR : public System
    {
      public:
        explicit SystemR(Parameters p, bool show_progress)
            : System(std::move(p)), m_show_progress(show_progress)
        {
        }

      protected:
        bool initRun() override
        {
            m_check_interval = std::max(2, parameters().sys.simulation_time / 100);
            m_last_percent   = -1;
            return true;
        }

        bool tearDownRun() override
        {
            if (m_show_progress) Rcpp::Rcout << '\n';
            return true;
        }

        void progressCallback(int t) override
        {
            if (!m_show_progress) return;
            const auto total = parameters().sys.simulation_time;
            if (total <= 0) return;
            const auto pct = std::min(100, (t * 100) / total);
            if (pct <= m_last_percent) return;
            if (pct == 100) {
                Rcpp::Rcout << "100%";
            } else if (pct % 10 == 0) {
                Rcpp::Rcout << pct << "%";
            } else {
                Rcpp::Rcout << '.';
            }
            m_last_percent = pct;
        }

        bool testForStop(int t) override
        {
            if (t % m_check_interval == 0)
            {
                Rcpp::checkUserInterrupt();
            }
            return false;
        }

      private:
        bool m_show_progress;
        int  m_check_interval = 2;
        int  m_last_percent   = -1;
    };

    Rcpp::List massSeriesToList(const std::vector<MassSeries>& comp_series,
                                const std::vector<std::string>& names,
                                const MassSeries& sink_series, const std::string& sink_name)
    {
        Rcpp::List out;
        for (std::size_t i = 0; i < comp_series.size(); ++i)
        {
            const auto& s = comp_series[i];
            if (!s.enabled) continue;
            Rcpp::List entry = Rcpp::List::create(Rcpp::Named("time")  = s.times,
                                                  Rcpp::Named("value") = s.values);
            out.push_back(entry, names[i]);
        }
        if (sink_series.enabled)
        {
            Rcpp::List entry = Rcpp::List::create(Rcpp::Named("time")  = sink_series.times,
                                                  Rcpp::Named("value") = sink_series.values);
            out.push_back(entry, sink_name);
        }
        return out;
    }

    Rcpp::List cdpToList(const std::vector<CdpSeries>& series,
                         const std::vector<std::string>& names)
    {
        Rcpp::List out;
        for (std::size_t i = 0; i < series.size(); ++i)
        {
            const auto& s = series[i];
            if (!s.enabled) continue;

            const auto n_t = s.times.size();
            const auto n_d = s.depths_um.size();

            Rcpp::NumericMatrix conc(static_cast<int>(n_d), static_cast<int>(n_t));
            for (std::size_t t = 0; t < n_t; ++t)
            {
                const auto& col = s.conc_per_time[t];
                for (std::size_t d = 0; d < n_d; ++d)
                {
                    conc(static_cast<int>(d), static_cast<int>(t)) = col[d];
                }
            }

            Rcpp::List entry = Rcpp::List::create(Rcpp::Named("time")     = s.times,
                                                  Rcpp::Named("depth_um") = s.depths_um,
                                                  Rcpp::Named("conc")     = conc);
            out.push_back(entry, names[i]);
        }
        return out;
    }

    Rcpp::List geometryToList(const Geometry& g)
    {
        return Rcpp::List::create(
            Rcpp::Named("min_step_um") = g.minSpaceStep(),
            Rcpp::Named("max_step_um") = g.maxSpaceStep(),
            Rcpp::Named("n_cells")     = g.size());
    }
}  // namespace

// [[Rcpp::export(name = ".cpp_validate", rng = false)]]
Rcpp::List cpp_validate(Rcpp::List params)
{
    Parameters p;
    try
    {
        p = parametersFromR(params);
    }
    catch (Rcpp::exception&) { throw; }
    catch (std::exception& e)
    {
        return Rcpp::List::create(Rcpp::Named("ok") = false,
                                  Rcpp::Named("error") = std::string(e.what()));
    }

    if (auto err = validate(p))
    {
        return Rcpp::List::create(Rcpp::Named("ok") = false, Rcpp::Named("error") = *err);
    }
    return Rcpp::List::create(Rcpp::Named("ok") = true,
                              Rcpp::Named("error") = R_NilValue);
}

// [[Rcpp::export(name = ".cpp_simulate", rng = false)]]
Rcpp::List cpp_simulate(Rcpp::List params, bool show_progress = false)
{
    Parameters p = parametersFromR(params);
    if (auto err = validate(p))
    {
        Rcpp::stop(*err);
    }

    SystemR sys(std::move(p), show_progress);
    const auto status = sys.run();

    std::string status_str = "executed";
    if (status == System::Result::Stopped) status_str = "stopped";
    if (status == System::Result::Failed)  status_str = "failed";

    const auto& parms = sys.parameters();
    return Rcpp::List::create(
        Rcpp::Named("status")   = status_str,
        Rcpp::Named("scaling")  = std::string(toString(parms.log.scaling)),
        Rcpp::Named("mass")     = massSeriesToList(sys.compartmentMass(),
                                                   sys.compartmentNames(),
                                                   sys.sinkMass(), parms.sink.name),
        Rcpp::Named("cdp")      = cdpToList(sys.cdp(), sys.compartmentNames()),
        Rcpp::Named("geometry") = geometryToList(sys.geometry()));
}
