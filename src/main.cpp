#include <mpi.h>

#include <iostream>

#include "muParser.h"
#include "parallelSolver.hpp"
#include "serialSolver.hpp"

using namespace laplace_solver;

int main(int argc, char* argv[]) {
    int ierr = MPI_Init(&argc, &argv);

    if (ierr != MPI_SUCCESS) {
        std::cerr << "Error initializing MPI! Error code: " << ierr
                  << std::endl;
        return 1;
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string filename;

    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "No input file name detected." << std::endl;
            std::cerr << "Using default filename \"data.json\"" << std::endl
                      << std::endl;
        }
        filename = "data.json";
    }
    else {
        filename = argv[1];
    }

    try {
        parallelSolver s(filename);
        s.solve();
        s.export_vtk("results.vtk");
    } catch (mu::Parser::exception_type& e) {
        std::cout << e.GetMsg() << std::endl;
    }

    ierr = MPI_Finalize();

    if (ierr != MPI_SUCCESS) {
        std::cerr << "Error finalizing MPI! Error code: " << ierr << std::endl;
        return 1;
    }

    return 0;
}
