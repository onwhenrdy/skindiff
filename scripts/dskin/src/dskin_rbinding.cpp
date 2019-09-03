#include "jsonparser.h"
#include "system.h"
#include <Rcpp.h>
#include <string>
// [[Rcpp::plugins("cpp11")]]

using namespace sc;

// R based System class
class SystemR : public System
{
  public:
    SystemR(const Parameter& parameter);

  protected:
    void progressCallback(int current_iteration) override final;
    bool testForStop(int current_iteration) override final;
    bool initRun() override final;
    bool tearDownRun() override final;

  private:
    int m_rcpp_check_interval;
    // progress counter
    int m_total_ticks;
    int m_last_perc;
};

SystemR::SystemR(const Parameter& parameter)
    : System(parameter), m_rcpp_check_interval(2), m_total_ticks(1), m_last_perc(-1)
{
}

void SystemR::progressCallback(int current_iteration)
{
    int percent = std::min(100, (current_iteration * 100) / m_total_ticks);
    if (percent <= m_last_perc)
    {
        return;
    }

    if (percent > m_last_perc)
    {
        if (percent >= 100)
        {
            Rcpp::Rcout << "100%\n";
        }
        else
        {
            Rcpp::Rcout << ((percent % 10 == 0) ? (std::to_string(percent) + "%") : ".");
        }
    }

    m_last_perc = percent;
}

bool SystemR::testForStop(int current_iteration)
{
    if (current_iteration % m_rcpp_check_interval == 0)
    {
        Rcpp::checkUserInterrupt();
    }

    return false;
}

bool SystemR::initRun()
{
    m_rcpp_check_interval = std::max(2, (this->simTime() / 100));
    return true;
}

bool SystemR::tearDownRun()
{
    m_last_perc = -1;
    return true;
}

///////////////////////////////////////////////

// helper functions
Rcpp::List writeToR(const Parameter& p, const SystemR& system, bool write_logs, bool success)
{
    auto result = Rcpp::List::create();

    if (p.logParameter().showProgressBar())
    {
        Rcpp::Rcout << "Build R data structure ...\n";
    }

    // meta data
    result.push_back(p.overviewString(), "log");
    result.push_back(p.logParameter().tag(), "tag");
    result.push_back(success, "success");

    auto geo_details = Rcpp::List::create();
    const auto& geo = system.geometry();
    geo_details.push_back(geo.minSpaceStep(), "min_ss");
    geo_details.push_back((geo.discMethod() == Geometry::DiscMethod::B_AND_K) ? geo.calculatedEta()
                                                                         : NA_REAL, "mb_eta");
    result.push_back(geo_details, "geometry");

    if (!write_logs)
    {
        return result;
    }

    // mass compartments
    const auto& layer_logger = system.compartmentLogger();
    const auto& s_logger     = system.sinkLogger();
    const auto& l_names      = p.layer_parameter();
    const auto v_name        = p.vehicleParameter().name();
    auto mass_df             = Rcpp::List::create();
    int i                    = 0;
    for (const auto& logger : layer_logger)
    {
        Rcpp::checkUserInterrupt();
        const std::string name = ((i == 0) ? v_name : l_names[i - 1].name());

        if (logger.enabled())
        {
            Rcpp::List tmp;
            tmp.push_back(logger.xs(), "time");
            tmp.push_back(logger.ys(), "mass");
            tmp.push_back(logger.filename(), "file.name");
            mass_df.push_back(tmp, name);
        }
        i++;
    }

    Rcpp::checkUserInterrupt();
    if (s_logger.enabled())
    {
        Rcpp::List tmp;
        tmp.push_back(s_logger.xs(), "time");
        tmp.push_back(s_logger.ys(), "mass");
        tmp.push_back(s_logger.filename(), "file.name");
        mass_df.push_back(tmp, p.sinkParameter().name());
    }
    mass_df.attr("class") = "data.frame";
    result.push_back(mass_df, "masses");

    // cdp data
    Rcpp::List conc_list;
    const auto& cdp_logger = system.cdpLogger();
    i                      = 0;
    for (const auto& logger : cdp_logger)
    {
        Rcpp::checkUserInterrupt();
        const std::string name = ((i == 0) ? v_name : l_names[i - 1].name());

        if (logger.enabled())
        {
            Rcpp::List tmp;

            const auto& d     = logger.data();
            const auto& times = logger.times();
            Rcpp::List data(d.size());
            Rcpp::CharacterVector data_names(times.size());
            int j = 0;
            for (const auto& column : d)
            {
                data_names[j] = std::to_string(times[j]);
                data[j]       = column;
                j++;
            }
            data.attr("names") = data_names;
            data.attr("class") = "data.frame";

            tmp.push_back(data, "data");
            tmp.push_back(logger.filename(), "file.name");
            tmp.push_back(times, "time");
            tmp.push_back(logger.space(), "x");
            conc_list.push_back(tmp, name);
        }
        i++;
    }

    result.push_back(conc_list, "concs");
    return result;
}

bool writeToFiles(const Parameter& p, const SystemR& system)
{
    bool ok = true;

    if (p.logParameter().showProgressBar())
    {
        Rcpp::Rcout << "Build logging files ...\n";
    }

    // compartment logs
    const auto s_logger = system.sinkLogger();
    if (s_logger.enabled())
    {
        ok = s_logger.writeToFile();
    }

    const auto layer_logger = system.compartmentLogger();
    if (ok)
    {
        for (const auto& logger : layer_logger)
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
    const auto cdp_logger = system.cdpLogger();
    if (ok)
    {
        for (const auto& logger : cdp_logger)
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


//' @keywords internal
// [[Rcpp::export(name = ".dskin.binding.geometry", rng = false)]]
Rcpp::List getGeometry(std::string json_parameter)
{
  auto result = Rcpp::List::create();

  JsonParser parser;
  auto ok = parser.parseFromString(json_parameter);
  if (!ok)
  {
    Rcpp::stop(parser.lastError());
    return result;
  }

  auto p = parser.parameter();
  SystemR system(p);

  const auto& geo = system.geometry();
  result.push_back(geo.minSpaceStep(), "min_ss");
  result.push_back((geo.discMethod() == Geometry::DiscMethod::B_AND_K) ? geo.calculatedEta()
                            : NA_REAL, "mb_eta");

  auto r_ss = Rcpp::List::create();

  const auto& space_steps = geo.spaceSteps();
  const auto& comps = system.compartments();
  for (const auto& c : comps)
  {
    const auto from_idx = c.geometryFromIdx();
    const auto to_idx = c.geometryToIdx();
    auto from = space_steps.begin() + from_idx;
    auto to = space_steps.begin() + to_idx;
    std::vector<double> local_ss(from, to + 1);
    r_ss.push_back(local_ss, c.name());
  }

  result.push_back(r_ss, "compartments");
  return result;
}


//' @keywords internal
// [[Rcpp::export(name = ".dskin.binding.simulate", rng = false)]]
Rcpp::List simulate(std::string json_parameter, bool write_to_R, bool write_to_files)
{
    JsonParser parser;
    auto ok = parser.parseFromString(json_parameter);
    if (!ok)
    {
        Rcpp::stop(parser.lastError());
    }
    else
    {
        auto p = parser.parameter();

        SystemR system(p);
        const auto result = system.run();
        if (result == System::Result::Executed)
        {
            if (write_to_files)
            {
                if (!writeToFiles(p, system))
                {
                    Rcpp::stop("Error: Could not write to log files.");
                }
            }

            return writeToR(p, system, write_to_R, true);
        }
        else
        {
            return writeToR(p, system, write_to_R, false);
        }
    }

    return Rcpp::List::create();
}
