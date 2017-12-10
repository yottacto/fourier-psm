#pragma once
#include "utils/type.hh"

namespace fpsm
{
    struct grid
    {
        grid(int ncd, int bnpcd);
        // return the 3d position
        point get_point(int rank, int local_id);
        // return the global 3d index
        index get_index(int rank, int local_id);

        // number of cores per direction
        int ncd{2};
        // base number of points per core per direction
        int bnpcd{3};
        // number of cores
        int nc{ncd * ncd * ncd};
        // number of points per core per direction
        int npcd{1 << bnpcd};
        // number of points per direction
        int npd{npcd * ncd};
        // number of points per core
        int npc{npcd * npcd * npcd};
        // number of points
        int np{npd * npd * npd};

    };
}

