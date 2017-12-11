#include <iostream>
#include "solver.hh"
#include "utils/constant.hh"

int main()
{

    // base number of cores per direction
    auto bncd = 1;
    // base number of points per core per direction
    auto bnpcd = 3;

    fpsm::solver s(bncd, bnpcd);
    s.init(fpsm::default_f);
    // s.iterate(100);
    // s.print();
}

