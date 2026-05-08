#include "geometry.h"
#include "parameter.h"
#include "system.h"

#include <testthat.h>

#include <cmath>
#include <numeric>
#include <vector>

using namespace sc;

namespace
{
    Parameters trivialParams(int sim_time = 30, int height_layer = 20)
    {
        Parameters p;
        p.sys.simulation_time = sim_time;
        p.sys.disc_method     = Geometry::DiscMethod::EQUI_DIST;
        p.sys.resolution      = 1;
        p.log.scaling         = Scaling::MG;
        p.vehicle.c_init   = 1.0;
        p.vehicle.height   = 30;
        p.vehicle.D        = 1.0;
        p.vehicle.app_area = 1.0;
        LayerParams sc;
        sc.name          = "SC";
        sc.height        = height_layer;
        sc.D             = 1.0;
        sc.K             = 1.0;
        sc.cross_section = 1.0;
        p.layers.push_back(sc);
        p.sink.Vd = 1.0;
        return p;
    }
}

context("Geometry")
{
    test_that("equidistant mesh assigns one cell per um per compartment")
    {
        Parameters p = trivialParams();
        std::vector<Compartment> comps;
        comps.push_back(Compartment{p.vehicle.height, p.vehicle.D, 1.0,
                                    p.vehicle.app_area * 1e8, p.vehicle.name});
        comps.push_back(Compartment{p.layers[0].height, p.layers[0].D, p.layers[0].K,
                                    p.vehicle.app_area * 1e8 * p.layers[0].cross_section,
                                    p.layers[0].name});
        Sink s;
        s.area_um2 = p.vehicle.app_area * 1e8;
        s.Vd       = 1.0;

        Geometry g;
        const auto ok = g.create(Geometry::DiscMethod::EQUI_DIST, comps, 1, &s);
        expect_true(ok);
        // 30 vehicle cells + 20 SC cells + 1 sink cell = 51
        expect_true(g.size() == 51);
        expect_true(comps[0].geo_from == 0);
        expect_true(comps[0].geo_to   == 29);
        expect_true(comps[1].geo_from == 30);
        expect_true(comps[1].geo_to   == 49);
        expect_true(s.geo_from == 50);
    }

    test_that("BK mesh adds refinement at compartment boundaries")
    {
        Parameters p = trivialParams();
        std::vector<Compartment> comps;
        comps.push_back(Compartment{p.vehicle.height, p.vehicle.D, 1.0,
                                    p.vehicle.app_area * 1e8, p.vehicle.name});
        comps.push_back(Compartment{p.layers[0].height, p.layers[0].D, p.layers[0].K,
                                    p.vehicle.app_area * 1e8 * p.layers[0].cross_section,
                                    p.layers[0].name});
        Sink s;
        s.area_um2 = p.vehicle.app_area * 1e8;

        Geometry g;
        g.setEta(0.6);
        const auto ok = g.create(Geometry::DiscMethod::B_AND_K, comps, 4, &s);
        expect_true(ok);
        expect_true(g.minSpaceStep() < 1.0);
        expect_true(g.maxSpaceStep() == 1.0);
    }
}

context("System mass conservation")
{
    test_that("perfect-sink runs conserve mass to ~machine epsilon")
    {
        Parameters p = trivialParams();
        p.vehicle.log_mass = true;
        p.layers[0].log_mass = true;
        p.sink.log_mass = true;
        System sys(std::move(p));
        const auto status = sys.run();
        expect_true(status == System::Result::Executed);

        const auto& m   = sys.compartmentMass();
        const auto& sm  = sys.sinkMass();
        expect_true(m.size() == 2);
        const auto n_pts = m[0].times.size();
        expect_true(n_pts == sm.times.size());

        const auto initial_total =
            m[0].values[0] + m[1].values[0] + sm.values[0];

        for (std::size_t i = 0; i < n_pts; ++i)
        {
            const auto total = m[0].values[i] + m[1].values[i] + sm.values[i];
            const auto rel   = std::abs(total - initial_total) / std::abs(initial_total);
            expect_true(rel < 1.0e-10);
        }
    }

    test_that("vehicle mass is monotonically non-increasing for K=1")
    {
        Parameters p = trivialParams();
        System sys(std::move(p));
        sys.run();
        const auto& v = sys.compartmentMass()[0].values;
        for (std::size_t i = 1; i < v.size(); ++i)
        {
            expect_true(v[i] <= v[i - 1] + 1e-12);
        }
    }

    test_that("sink mass is monotonically non-decreasing for a perfect sink")
    {
        Parameters p = trivialParams();
        System sys(std::move(p));
        sys.run();
        const auto& s = sys.sinkMass().values;
        for (std::size_t i = 1; i < s.size(); ++i)
        {
            expect_true(s[i] + 1e-12 >= s[i - 1]);
        }
    }

    test_that("CDP profiles are recorded when enabled")
    {
        Parameters p = trivialParams();
        p.vehicle.log_cdp = true;
        p.layers[0].log_cdp = true;
        System sys(std::move(p));
        sys.run();

        const auto& cdp = sys.cdp();
        expect_true(cdp.size() == 2);
        for (const auto& s : cdp)
        {
            expect_true(s.enabled);
            expect_false(s.depths_um.empty());
            expect_false(s.times.empty());
            expect_true(s.times.size() == s.conc_per_time.size());
            for (const auto& profile : s.conc_per_time)
            {
                expect_true(profile.size() == s.depths_um.size());
            }
        }
    }
}

context("Parameter validation")
{
    test_that("default Parameters is valid (no layers, single vehicle)")
    {
        Parameters p;
        // Default vehicle.height = 10, c_init = 1, app_area = 1, D = 1 is valid.
        // No layers, vehicle.removed() is false -> valid.
        expect_false(static_cast<bool>(validate(p)));
    }

    test_that("removing the vehicle without any layer is rejected")
    {
        Parameters p;
        p.vehicle.remove_at = 60;
        const auto err = validate(p);
        expect_true(static_cast<bool>(err));
    }

    test_that("layer with cross_section > 1 is rejected")
    {
        Parameters p;
        LayerParams l;
        l.name = "bad";
        l.cross_section = 2.0;
        p.layers.push_back(l);
        const auto err = validate(p);
        expect_true(static_cast<bool>(err));
    }
}
