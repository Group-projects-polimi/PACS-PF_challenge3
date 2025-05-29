#ifndef SERIAL_SOLVER_HPP
#define SERIAL_SOLVER_HPP

#include <memory>
#include <vector>

#include "matrix.hpp"
#include "parameters.hpp"

using namespace algebra;

namespace laplace_solver {

/**
 * @class serialSolver
 * @brief Serial solver for the Laplacian equation on a square domain using the Jacobi method.
 *
 * This class manages the serial solution of the Laplace equation, including mesh generation,
 * parameter handling, solution initialization, and result export.
 */
class serialSolver {
   private:
    Parameters param;                ///< Problem parameters.
    std::unique_ptr<matrix> solution;///< Solution matrix.
    std::vector<double> mesh_x;      ///< X-coordinates of the mesh.
    std::vector<double> mesh_y;      ///< Y-coordinates of the mesh.
    double mesh_size;                ///< Grid spacing.

   public:
    /**
     * @brief Constructs the solver with given parameters.
     * @param p Parameters for the Laplacian problem.
     */
    serialSolver(Parameters const& p) : param{p} {
        create_mesh();
        initialize_sol();
    }

    /**
     * @brief Constructs the solver by reading parameters from a file.
     * @param filename Path to the parameter file.
     */
    serialSolver(std::string const& filename) {
        readParameters(filename);
        create_mesh();
        initialize_sol();
    }

    /**
     * @brief Constructs the solver by reading parameters from a file and setting the number of grid points.
     * @param filename Path to the parameter file.
     * @param n Number of grid points per spatial direction.
     */
    serialSolver(std::string const& filename, size_t n) {
        readParameters(filename);
        serial_setpts(n);
        create_mesh();
        initialize_sol();
    }

    /**
     * @brief Reads problem parameters from a file.
     * @param filename Path to the parameter file.
     */
    void readParameters(std::string const& filename);

    /**
     * @brief Creates the computational mesh for the domain.
     */
    void create_mesh();

    /**
     * @brief Initializes the solution matrix.
     */
    void initialize_sol();

    /**
     * @brief Runs the Jacobi solver.
     */
    void solve();

    /**
     * @brief Prints the error between the computed and exact solution (if available).
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
     * @brief Sets the number of grid points per direction.
     * @param n Number of grid points.
     */
    void serial_setpts(size_t n);
};

}  // namespace laplace_solver

#endif
