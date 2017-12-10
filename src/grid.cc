#pragma once
#include "grid.hh"

namespace fpsm
{
    grid::grid(int bncd, int bnpcd) :
        bncd{bncd},
        ncd{1 << bncd},
        bnpcd{bnpcd},
        nc{ncd * ncd * ncd},
        npcd{1 << bnpcd},
        npd{npcd * ncd},
        npc{npcd * npcd * npcd},
        np{npd * npd * npd},
        delta{len / npd}
    {
    }

    index grid::get_core_index(int rank) const
    {
        auto z  = rank / (ncd * ncd);
        auto y  = (rank % (ncd * ncd)) / ncd;
        auto x  = (rank % (ncd * ncd)) % ncd;
        return {x, y, z};
    }

    index grid::get_index(int rank, int local_id) const
    {
        auto cindex = get_core_index(rank);

        auto lz = local_id / (npcd * npcd);
        auto ly = (local_id % (npcd * npcd)) / npcd;
        auto lx = (local_id % (npcd * npcd)) % npcd;

        return {
            cindex.x * npcd + lx,
            cindex.y * npcd + ly,
            cindex.z * npcd + lz
        };
    }

    point grid::get_point(int rank, int local_id) const
    {
        auto index = get_index(rank, local_id);
        return {
            -2. * pi + index.x * delta,
            -2. * pi + index.y * delta,
            -2. * pi + index.z * delta
        };
    }

    int grid::get_core_rank(index const& p) const
    {
        return p.z * ncd * ncd + p.y * ncd + p.z;
    }

    int grid::get_local_id(index const& p) const
    {
        return p.z * npcd * npcd + p.y * npcd + p.z;
    }
}

