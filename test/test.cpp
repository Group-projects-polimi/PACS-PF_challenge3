#include "test.hpp"

#include <mpi.h>

#include <chrono>
#include <parallelSolver.hpp>
#include <serialSolver.hpp>
#include <string>

using namespace laplace_solver;

/**
 * @brief Runs serial performance tests for the Jacobi Laplacian solver.
 * @details
 * For grid sizes from 5 to 100 (step 5), constructs a serial solver, times the solution,
 * and writes the grid size and elapsed time (in milliseconds) to a file.
 * The output file is "./test/output/serialtime.txt".
 * @param filename Path to the parameter file.
 */
void run_serialtest(std::string filename) {
    std::ofstream output_serial{"./test/output/serialtime.txt"};

    for (size_t n = 10; n <= 100; n += 5) {
        serialSolver s1(filename, n);
        std::chrono::duration<double, std::milli> duration1;
        auto start = std::chrono::high_resolution_clock::now();
        s1.solve();
        auto end = std::chrono::high_resolution_clock::now();
        duration1 = end - start;
        output_serial << n << " " << duration1.count() << std::endl;
    }
}

/**
 * @brief Runs parallel performance tests for the Jacobi Laplacian solver using MPI.
 * @details
 * For grid sizes from 5 to 100 (step 5), constructs a parallel solver, times the solution
 * (only on rank 0), and writes the grid size and elapsed time (in milliseconds) to a file.
 * The output file is "./test/output/paralleltime_<size>.txt", where <size> is the number of MPI ranks.
 * @param filename Path to the parameter file.
 */
void run_paralleltest(std::string filename) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::ofstream output_parallel;
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double, std::milli> duration1;

    if (rank == 0) {
        output_parallel.open("./test/output/paralleltime_" +
                             std::to_string(size) + ".txt");
    }

    for (size_t n = 10; n <= 100; n += 5) {
        parallelSolver s1(filename, n);
        if (rank == 0) {
            start = std::chrono::high_resolution_clock::now();
        }
        s1.solve();
        if (rank == 0) {
            end = std::chrono::high_resolution_clock::now();
            duration1 = end - start;
            output_parallel << n << " " << duration1.count() << std::endl;
        }
    }
}
