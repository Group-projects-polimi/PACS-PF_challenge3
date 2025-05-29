#include <iostream>

#include "muParser.h"
#include "parallelSolver.hpp"
#include "serialSolver.hpp"
#include "test.hpp"

using namespace laplace_solver;

#ifdef PARALLEL
#include <mpi.h>
int main(int argc, char* argv[]) {
    int ierr = MPI_Init(&argc, &argv);

    if (ierr != MPI_SUCCESS) {
        std::cerr << "Error initializing MPI! Error code: " << ierr
                  << std::endl;
        return 1;
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string filename = argv[1];
    try {
        run_paralleltest(filename);
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
#else
int main(int argc, char* argv[]) {
    std::string filename = argv[1];

    try {
        run_serialtest(filename);
    } catch (mu::Parser::exception_type& e) {
        std::cout << e.GetMsg() << std::endl;
    }

    return 0;
}
#endif
