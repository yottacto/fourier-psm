#include <iostream>
#include "solver.hh"

int main()
{

    // number of cores per direction
    auto ncd = 2;
    // base number of points per core per direction
    auto bnpcd = 3;
    // number of cores
    auto nc = ncd * ncd * ncd;
    // number of points per core per direction
    auto npcd = 1 << bnpcd;
    // number of points per direction
    auto npd = npcd * ncd;
    // number of points per core
    auto npc = npcd * npcd * npcd;
    // number of points
    auto np = npd * npd * npd;

    fpsm::solver s(ncd, bnpcd);
    s.init();
    s.iterate(100);
    s.print();
}

