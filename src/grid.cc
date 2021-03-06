#include <iostream>
#include "grid.hh"

namespace fpsm
{
    grid::grid(int bncd, int bnpcd) :
        bncd{bncd},
        ncd{1 << bncd},
        nc{ncd * ncd * ncd},
        bnpcd{bnpcd},
        npcd{1 << bnpcd},
        npd{npcd * ncd},
        npc{npcd * npcd * npcd},
        np{npd * npd * npd},
        delta{len / npd}
    {
    }

    void grid::print() const
    {
        std::cerr << "grid:\n";
        std::cerr << "bncd = " << bncd << "\n";
        std::cerr << "ncd = " << ncd << "\n";
        std::cerr << "nc = " << nc << "\n";
        std::cerr << "bnpcd = " << bnpcd << "\n";
        std::cerr << "npcd = " << npcd << "\n";
        std::cerr << "npd = " << npd << "\n";
        std::cerr << "npc = " << npc << "\n";
        std::cerr << "np = " << np << "\n";
    }

    index grid::get_index(int global_id) const
    {
        auto z = global_id / (npd * npd);
        auto y = (global_id % (npd * npd)) / npd;
        auto x = (global_id % (npd * npd)) % npd;
        return {x, y, z};
    }

    index grid::get_core_index(int rank) const
    {
        auto z  = rank / (ncd * ncd);
        auto y  = (rank % (ncd * ncd)) / ncd;
        auto x  = (rank % (ncd * ncd)) % ncd;
        return {x, y, z};
    }

    index grid::get_local_index(int local_id) const
    {
        auto z  = local_id / (npcd * npcd);
        auto y  = (local_id % (npcd * npcd)) / npcd;
        auto x  = (local_id % (npcd * npcd)) % npcd;
        return {x, y, z};
    }

    index grid::get_index(int rank, int local_id) const
    {
        auto cindex = get_core_index(rank);
        auto lindex = get_local_index(local_id);

        return {
            cindex.x * npcd + lindex.x,
            cindex.y * npcd + lindex.y,
            cindex.z * npcd + lindex.z
        };
    }

    point grid::get_point(int global_id) const
    {
        auto index = get_index(global_id);
        return {
            -pi + index.x * delta,
            -pi + index.y * delta,
            -pi + index.z * delta
        };
    }

    point grid::get_point(int rank, int local_id) const
    {
        auto index = get_index(rank, local_id);
        return {
            -pi + index.x * delta,
            -pi + index.y * delta,
            -pi + index.z * delta
        };
    }

    int grid::get_core_rank(index const& p) const
    {
        return p.z * ncd * ncd + p.y * ncd + p.x;
    }

    int grid::get_local_id(index const& p) const
    {
        return p.z * npcd * npcd + p.y * npcd + p.x;
    }

    int grid::get_id(index const& p) const
    {
        return p.z * npd * npd + p.y * npd + p.x;
    }

    int grid::get_id(int rank, int local_id) const
    {
        auto index = get_index(rank, local_id);
        return get_id(index);
    }
}

