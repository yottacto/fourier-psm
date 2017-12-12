#include <iostream>
#include <mpi.h>
#include "solver.hh"
#include "utils/constant.hh"

#include "fft/fft.hh"
#include "utils/type.hh"

int main()
{
    MPI::Init();
    // base number of cores per direction
    auto bncd = 2;
    // base number of points per core per direction
    auto bnpcd = 3;

    fpsm::solver s(bncd, bnpcd);

    s.init(fpsm::debug_f);

    // s.init(fpsm::default_f);
    // s.iterate(1);
    s.print();

    MPI::Finalize();
}

