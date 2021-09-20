/**********************************************************************
 *                     Unit testing for Protein code                      *
 **********************************************************************/

#include "ProteinBindStatus.hpp"
#include "ProteinConfig.hpp"
#include "ProteinData.hpp"
#include "ProteinType.hpp"

// SimToolbox submodule
#include "FDPS/particle_simulator.hpp"
#include "Sylinder/SylinderNear.hpp"

#include "SRC/catch.hpp"

#include <array>
#include <cassert>

SylinderNearEP MockSylinder(int id) {
    SylinderNearEP rod;
    rod.gid = id;
    rod.length = 400;
    rod.rank = 0;
    rod.radius = .5;
    for (int i = 0; i < 3; ++i) {
        rod.pos[i] = 0;
        rod.direction[i] = 0;
    }
    rod.direction[0] = 1;
    return rod;
}

TEST_CASE("Parallel to Anti-parallel binding affinity", "[PtoAPratioSD]") {
    constexpr double small = 1e-6;
    ProteinData pData;
    pData.setMockProtein();
    double PtoAPratio = .5;
    pData.property.PtoAPratio = PtoAPratio;
    SylinderNearEP rod = MockSylinder(0);
    double dt = .001;
    int end_bound = 0;
    double end_loc = 0;
    pData.bind.setBind(end_bound, rod.gid, rod.globalIndex, rod.direction,
                       rod.pos, end_loc, rod.length, rod.rank);

    const double D = 0.024;
    const double alpha = 0.1 / (2 * 0.00411);
    const double freelength = 0.05;
    const double M = alpha * D * D;
    const double ell0 = freelength / D;

    double grid_num = 1024;
    LUTFillerEdep *lut_filler_ptr = new LUTFillerEdep(grid_num, grid_num);
    lut_filler_ptr->Init(M, ell0, D);

    LookupTable LUT(lut_filler_ptr);

    pData.property.LUTablePtr = &LUT;

    SECTION("Check binding ratio") {
        double paraBindFactor =
            pData.getBindingFactorSD(1 - end_bound, rod.direction);
        double antiparaDir[3];
        for (int i = 0; i < 3; ++i) {
            antiparaDir[i] = -1. * rod.direction[i];
        }
        double antiparaBindFactor =
            pData.getBindingFactorSD(1 - end_bound, antiparaDir);
        CHECK(abs(PtoAPratio - (paraBindFactor / antiparaBindFactor)) < small);
    }
}
