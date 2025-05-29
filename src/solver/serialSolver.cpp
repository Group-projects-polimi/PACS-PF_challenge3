#include "serialSolver.hpp"

#include <mpi.h>

#include <fstream>
#include <iostream>
#include <string>

#include "json.hpp"
#include "muParser.h"

namespace laplace_solver {

/**
 * @brief Creates the computational mesh for the domain.
 * @details
 * Initializes the mesh_x and mesh_y vectors with the coordinates of the grid points
 * in the x and y directions, respectively. The mesh is uniform and spans [0,1] in both directions.
 */
void serialSolver::create_mesh() {
    mesh_x.resize(param.grid_pts);
    mesh_y.resize(param.grid_pts);

    mesh_x[0] = 0;
    mesh_y[0] = 0;
    mesh_x[param.grid_pts - 1] = 1.0;
    mesh_y[param.grid_pts - 1] = 1.0;

    mesh_size = 1.0 / (param.grid_pts - 1);

    for (size_t i = 1; i < param.grid_pts - 1; ++i) {
        mesh_x[i] = mesh_x[i - 1] + mesh_size;
        mesh_y[i] = mesh_y[i - 1] + mesh_size;
    }
}

/**
 * @brief Initializes the solution matrix with boundary conditions.
 * @details
 * Allocates the solution matrix and sets the boundary values using the muParser
 * boundary condition. Handles all four boundaries (top, bottom, left, right) and corners.
 */
void serialSolver::initialize_sol() {
    solution = std::make_unique<matrix>(param.grid_pts, param.grid_pts);

    // Set the boundary conditions at the 4 corners of the mesh outside the
    // loop so that each corner is evaluated only once
    param.boundary_point[0] = mesh_x[0];
    param.boundary_point[1] = mesh_y[0];
    (*solution)(param.grid_pts - 1, 0) = param.boundary_condition.Eval();

    param.boundary_point[0] = mesh_x[param.grid_pts - 1];
    param.boundary_point[1] = mesh_y[0];
    (*solution)(param.grid_pts - 1, param.grid_pts - 1) =
        param.boundary_condition.Eval();

    param.boundary_point[0] = mesh_x[0];
    param.boundary_point[1] = mesh_y[param.grid_pts - 1];
    (*solution)(0, 0) = param.boundary_condition.Eval();

    param.boundary_point[0] = mesh_x[param.grid_pts - 1];
    param.boundary_point[1] = mesh_y[param.grid_pts - 1];
    (*solution)(0, param.grid_pts - 1) = param.boundary_condition.Eval();

    for (size_t i = 1; i < param.grid_pts - 1; ++i) {
        param.boundary_point[0] = mesh_x[i];
        param.boundary_point[1] = mesh_y[0];
        (*solution)(param.grid_pts - 1, i) = param.boundary_condition.Eval();

        param.boundary_point[0] = mesh_x[i];
        param.boundary_point[1] = mesh_y[param.grid_pts - 1];
        (*solution)(0, i) = param.boundary_condition.Eval();

        param.boundary_point[0] = mesh_x[0];
        param.boundary_point[1] = mesh_y[i];
        (*solution)(param.grid_pts - 1 - i, 0) =
            param.boundary_condition.Eval();

        param.boundary_point[0] = mesh_x[param.grid_pts - 1];
        param.boundary_point[1] = mesh_y[i];
        (*solution)(param.grid_pts - 1 - i, param.grid_pts - 1) =
            param.boundary_condition.Eval();
    }
}

/**
 * @brief Prints the current problem parameters.
 * @details
 * Uses the overloaded operator<< to print all parameters to the standard output.
 */
void serialSolver::print_parameters() const { std::cout << param << std::endl; }

/**
 * @brief Prints the mesh coordinates.
 * @details
 * Prints the (x, y) coordinates of each grid node in a formatted table.
 */
void serialSolver::print_mesh() const {
    std::cout << "PRINTING THE MESH" << std::endl;
    for (long long i = param.grid_pts - 1; i >= 0; --i) {
        for (long long j = 0; j < param.grid_pts; ++j) {
            std::stringstream ss;
            ss << "(" << mesh_x[j] << ", " << mesh_y[i] << ")";
            std::string s = ss.str();
            std::cout << std::setw(14) << s;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/**
 * @brief Prints the computed solution.
 * @details
 * Prints the solution matrix in a formatted table.
 */
void serialSolver::print_solution() const {
    std::cout << "PRINTING THE SOLUTION" << std::endl;
    for (long long i = 0; i < param.grid_pts; ++i) {
        for (long long j = 0; j < param.grid_pts; ++j) {
            std::stringstream ss;
            ss << (*solution)(i, j);
            std::string s = ss.str();
            std::cout << std::setw(14) << s;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/**
 * @brief Reads problem parameters from a file.
 * @details
 * Reads a JSON file containing all parameters and expressions. If the file is not found,
 * assigns default values. Sets up muParser objects and defines variables.
 * @param filename Path to the parameter file.
 */
void serialSolver::readParameters(std::string const& filename) {
    std::ifstream jfile(filename);

    if (!jfile) {
        std::cerr << "ERROR: parameters file " << filename << " does not exist."
                  << std::endl;
        std::cerr << "Reverting to default values." << std::endl;
        param.assign_defaults();
        return;
    }

    nlohmann::json data = nlohmann::json::parse(jfile);
    param.max_it = data["max_iterations"].get<unsigned long long>();
    param.tolerance = data["tolerance"].get<double>();
    param.grid_pts =
        data["grid_points_in_each_direction"].get<unsigned long long>();
    param.exact_sol_given = data["exact_solution_given"].get<bool>();
    std::string fun = data["forcing_term"].get<std::string>();
    std::string bc = data["boundary_condition"].get<std::string>();
    std::string exact = data["exact_solution"].get<std::string>();
    jfile.close();

    mu::Parser forcing_term, boundary_condition, exact_solution;
    forcing_term.SetExpr(fun);
    param.forcing_term = forcing_term;

    boundary_condition.SetExpr(bc);
    param.boundary_condition = boundary_condition;

    exact_solution.SetExpr(exact);
    param.exact_solution = exact_solution;

    param.define_variables();
}

/**
 * @brief Runs the Jacobi solver.
 * @details
 * Iteratively updates the solution matrix using the Jacobi method until the convergence
 * criterion is met or the maximum number of iterations is reached. The force term is
 * precomputed for efficiency.
 */
void serialSolver::solve() {
    std::unique_ptr<matrix> new_solution =
        std::make_unique<matrix>(param.grid_pts, param.grid_pts);

    // Copy boundary conditions outside the solver loop so
    // that it's done only once
    (*new_solution)(0, 0) = (*solution)(0, 0);
    (*new_solution)(param.grid_pts - 1, 0) = (*solution)(param.grid_pts - 1, 0);
    (*new_solution)(0, param.grid_pts - 1) = (*solution)(0, param.grid_pts - 1);
    (*new_solution)(param.grid_pts - 1, param.grid_pts - 1) =
        (*solution)(param.grid_pts - 1, param.grid_pts - 1);

    for (size_t i = 1; i < param.grid_pts - 1; ++i) {
        (*new_solution)(i, 0) = (*solution)(i, 0);
        (*new_solution)(i, param.grid_pts - 1) =
            (*solution)(i, param.grid_pts - 1);
    }

    for (size_t i = 1; i < param.grid_pts - 1; ++i) {
        (*new_solution)(0, i) = (*solution)(0, i);
        (*new_solution)(param.grid_pts - 1, i) =
            (*solution)(param.grid_pts - 1, i);
    }

    // Build the force term
    matrix force(param.grid_pts - 2, param.grid_pts - 2);
    for (size_t i = 0; i < param.grid_pts - 2; ++i) {
        for (size_t j = 0; j < param.grid_pts - 2; ++j) {
            param.force_point[0] = mesh_x[j + 1];
            param.force_point[1] = mesh_y[i + 1];
            force(param.grid_pts - 3 - i, j) =
                mesh_size * mesh_size * param.forcing_term.Eval();
        }
    }

    size_t k = 0;
    double conv_crit = param.tolerance + 1.0;
    while (k < param.max_it && conv_crit > param.tolerance) {
        conv_crit = 0.0;
        for (size_t j = 1; j < param.grid_pts - 1; ++j) {
            for (size_t i = 1; i < param.grid_pts - 1; ++i) {
                (*new_solution)(i, j) =
                    0.25 * ((*solution)(i - 1, j) + (*solution)(i + 1, j) +
                            (*solution)(i, j - 1) + (*solution)(i, j + 1) +
                            force(i - 1, j - 1));
                conv_crit += ((*new_solution)(i, j) - (*solution)(i, j)) *
                             ((*new_solution)(i, j) - (*solution)(i, j));
            }
        }
        conv_crit = std::sqrt(mesh_size * conv_crit);
        std::swap(solution, new_solution);
        ++k;
    }

#ifdef DEBUG
    if (k == param.max_it)
        std::cout << "Warning: serial solver reached max number of iterations "
                     "without "
                     "satisfying the convergence criterion, stopping."
                  << std::endl;
#endif
}

/**
 * @brief Computes and prints the error between the computed and exact solution.
 * @details
 * For each grid node, evaluates the exact solution and computes the squared difference
 * with the computed solution. The square root of the sum is printed.
 */
void serialSolver::print_error() {
    double error = 0.0;

    for (size_t i = 1; i < param.grid_pts - 1; ++i) {
        for (size_t j = 1; j < param.grid_pts - 1; ++j) {
            param.exact_point[0] = mesh_x[j];
            param.exact_point[1] = mesh_y[i];
            double exact = param.exact_solution.Eval();
            error += (exact - (*solution)(j, param.grid_pts - 1 - i)) *
                     (exact - (*solution)(j, param.grid_pts - 1 - i));
        }
    }

    std::cout << "The square root of the sum of the squared differences "
                 "between the computeted solution and the provided exact "
                 "solution in the grid nodes is: "
              << std::sqrt(error) << std::endl
              << std::endl;
}

/**
 * @brief Sets the number of grid points per direction.
 * @param n Number of grid points.
 */
void serialSolver::serial_setpts(size_t n) { param.grid_pts = n; }

/**
 * @brief Exports the solution to a VTK file for visualization.
 * @details
 * Writes the solution matrix to a VTK file in ASCII format for visualization.
 * @param filename Output VTK file name.
 * @return True if export was successful, false otherwise.
 */
bool serialSolver::export_vtk(std::string filename) const {
    std::ofstream vtkfile{filename};

    vtkfile << "# vtk DataFile Version 3.0" << std::endl;
    vtkfile << "Laplace Solver Solution" << std::endl;
    vtkfile << "ASCII" << std::endl;
    vtkfile << "DATASET STRUCTURED_POINTS" << std::endl;
    vtkfile << "DIMENSIONS " << param.grid_pts << " " << param.grid_pts << " "
            << "1" << std::endl;
    vtkfile << "ORIGIN 0.0 0.0 0.0" << std::endl;
    vtkfile << "SPACING " << mesh_size << " " << mesh_size << " " << "1.0"
            << std::endl;
    vtkfile << "POINT_DATA " << param.grid_pts * param.grid_pts << std::endl;
    vtkfile << "SCALARS solution_U double 1" << std::endl;
    vtkfile << "LOOKUP_TABLE default" << std::endl;

    for (int j = 0; j < param.grid_pts; ++j) {
        for (int i = 0; i < param.grid_pts; ++i) {
            vtkfile << (*solution)(j, i) << std::endl;
        }
    }
    return true;
}

}  // namespace laplace_solver
