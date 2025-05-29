#ifndef PARALLEL_SOLVER_HPP
#define PARALLEL_SOLVER_HPP

#include <mpi.h>

#include <memory>
#include <vector>

#include "matrix.hpp"
#include "parameters.hpp"

using namespace algebra;

namespace laplace_solver {

/**
 * @class parallelSolver
 * @brief Parallel solver for the Laplacian equation on a square domain using the Jacobi method and MPI/OpenMP.
 *
 * This class manages the parallel solution of the Laplace equation, including mesh generation,
 * parameter handling, solution initialization, and result export. Each MPI process handles a portion
 * of the domain.
 */
class parallelSolver {
   private:
    int rank;        ///< MPI rank of the current process.
    int size;        ///< Total number of MPI processes.
    int local_n;     ///< Number of rows handled by this process.
    int global_row;  ///< Global row index offset for this process.
    bool last_rank;  ///< True if this process is the last in rank order.
    bool first_rank; ///< True if this process is the first in rank order.
    double mesh_size;///< Grid spacing.

    std::unique_ptr<matrix> solution; ///< Local solution matrix.
    std::vector<double> mesh_x;       ///< X-coordinates of the mesh.
    std::vector<double> mesh_y;       ///< Y-coordinates of the mesh.
    Parameters param;                 ///< Problem parameters.

   public:
    /**
     * @brief Constructs the solver with given parameters.
     * @param p Parameters for the Laplacian problem.
     */
    parallelSolver(Parameters const& p) : param{p} {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        initialize();
        create_mesh();
        initialize_sol();
    }

    /**
     * @brief Constructs the solver by reading parameters from a file.
     * @param filename Path to the parameter file.
     */
    parallelSolver(std::string const& filename) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        readParameters(filename);
        initialize();
        create_mesh();
        initialize_sol();
    }

    /**
     * @brief Constructs the solver by reading parameters from a file and setting the number of grid points.
     * @param filename Path to the parameter file.
     * @param n Number of grid points per spatial direction.
     */
    parallelSolver(std::string const& filename, size_t n) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        readParameters(filename);
        parallel_setpts(n);
        initialize();
        create_mesh();
        initialize_sol();
    }

    /**
     * @brief Initializes MPI-related and solver-specific variables.
     */
    void initialize();

    /**
     * @brief Reads parameters from a file and initializes the solver.
     * @param filename Path to the parameter file.
     */
    void readParameters(std::string const& filename);

    /**
     * @brief Creates the computational mesh for the domain.
     */
    void create_mesh();

    /**
     * @brief Initializes the solution matrix with boundary conditions and forcing terms.
     */
    void initialize_sol();

    /**
     * @brief Solves the Laplace equation using the Jacobi method.
     */
    void solve();

    /**
     * @brief Computes the error between the computed solution and the exact solution (if provided).
     */
    void print_error();

    /**
     * @brief Prints the current problem parameters.
     */
    void print_parameters() const;

    /**
     * @brief Prints the mesh coordinates.
     */
    void print_mesh() const;

    /**
     * @brief Prints the computed solution.
     */
    void print_solution() const;

    /**
     * @brief Exports the solution to a VTK file for visualization.
     * @param filename Output VTK file name.
     * @return True if export was successful, false otherwise.
     */
    bool export_vtk(std::string filename) const;

    /**
     * @brief Sets the number of grid points per direction in parallel.
     * @param n Number of grid points.
     */
    void parallel_setpts(size_t n);
};

}  // namespace laplace_solver

#endif
