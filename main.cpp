#include "helper.h"
#include "macros.h"
#include "session.h"
#include <iostream>

using namespace sc;

void bench_test(const VersionInfo& v_info, int times, int nargs, const char* args[]);

void bench_test(const VersionInfo& v_info, int times, int nargs, const char* args[])
{
    for (int i = 0; i < times; ++i)
    {
        Session session(v_info);
        session.runFromCmdLineArgs(nargs, args);
    }
}

int main(int nargs, const char* args[])
{
    VersionInfo v_info(STRING(DSCMDAPPNAME), DSCMDMAJVERSION, DSCMDMINVERSION, DSCMDPATCHLEVEL);
    v_info.setBuildId(STRING(DSCMDBUILDNUMBER));
    v_info.setFullName("The DSkin Command Line Tool");
    v_info.setCopyrightNote("Scientific Consilience GmbH");

    /*
    using MESURE = bench::measure<std::chrono::milliseconds, std::chrono::system_clock>;
    static const auto N_RUNS = 50;
    auto run_time = MESURE::duration(bench_test, v_info, N_RUNS, nargs, args);
    std::cerr << "Duration [ms]: " << static_cast<double>(run_time.count())/N_RUNS;
    auto ok = true;
    */

    Session session(v_info);
    auto ok = session.runFromCmdLineArgs(nargs, args);
    return (ok) ? 0 : 1;
}
