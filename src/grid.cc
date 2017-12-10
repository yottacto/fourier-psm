#pragma once
#include "grid.hh"

namespace fpsm
{
    grid::grid(int ncd, int bnpcd) :
        ncd{ncd},
        bnpcd{bnpcd},
        nc{ncd * ncd * ncd},
        npcd{1 << bnpcd},
        npd{npcd * ncd},
        npc{npcd * npcd * npcd},
        np{npd * npd * npd}
    {
    }
}

