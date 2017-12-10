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
        np{npd * npd * npd}
    {
    }
}

