#pragma once
#include "utils/type.hh"
#include "utils/constant.hh"

namespace fpsm
{
    struct grid
    {
        grid() = default;
        grid(int bncd, int bnpcd);

        // return the 3d position
        point get_point(int rank, int local_id) const;

        // return the 3d position
        point get_point(int global_id) const;

        // return the global 3d index
        index get_index(int rank, int local_id) const;

        // return the global 3d index
        index get_index(int global_id) const;

        // return the global 3d index of cores
        index get_core_index(int rank) const;

        // return the local 3d index of cores
        index get_local_index(int local_id) const;

        int get_core_rank(index const& p) const;
        int get_local_id(index const& p) const;
        int get_id(index const& p) const;
        int get_id(int rank, int local_id) const;

        void print() const;


        // base number of cores per direction
        int bncd{1};

        // number of cores per direction
        int ncd{1 << bncd};

        // number of cores
        int nc{ncd * ncd * ncd};

        // base number of points per core per direction
        int bnpcd{3};

        // number of points per core per direction
        int npcd{1 << bnpcd};

        // number of points per direction
        int npd{npcd * ncd};

        // number of points per core
        int npc{npcd * npcd * npcd};

        // number of points
        int np{npd * npd * npd};

        double len{2. * pi};
        double delta{len / npd};
    };
}

