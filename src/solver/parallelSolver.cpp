#include "parallelSolver.hpp"

#include <omp.h>

#include <format>
#include <fstream>
#include <iostream>
#include <string>

#include "json.hpp"
#include "muParser.h"

namespace laplace_solver {

/**
 * @brief Initializes MPI-related and solver-specific variables.
 * @details
 * Computes the number of rows handled by this process (local_n), the global row offset (global_row),
 * and sets flags for first and last rank. Uses MPI_Exscan to determine the offset for each process.
 */
void parallelSolver::initialize() {
    local_n = (param.grid_pts / size) + (rank < (param.grid_pts % size));

    global_row = 0;
    MPI_Exscan(&local_n, &global_row, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    first_rank = (rank == 0);
    last_rank = (rank == size - 1);
}

/**
 * @brief Creates the computational mesh for the domain.
 * @details
 * Initializes the mesh_x and mesh_y vectors with the coordinates of the grid points
 * in the x and y directions, respectively. Handles the assignment of boundary points
 * for the first and last ranks.
 */
void parallelSolver::create_mesh() {
    mesh_x.resize(param.grid_pts);
    mesh_y.resize(local_n);

    mesh_size = 1.0 / (param.grid_pts - 1);

    mesh_x[0] = 0;
    mesh_y[0] = mesh_size * global_row;
    mesh_x[param.grid_pts - 1] = 1.0;

    if (last_rank) mesh_y[local_n - 1] = 1.0;

    for (size_t i = 1; i < param.grid_pts - 1; ++i) {
        mesh_x[i] = mesh_x[i - 1] + mesh_size;
    }

    for (size_t i = 1; i < local_n - last_rank; ++i) {
        mesh_y[i] = mesh_y[i - 1] + mesh_size;
    }
}

/**
 * @brief Initializes the solution matrix with boundary conditions.
 * @details
 * Allocates the local solution matrix and sets the boundary values using the muParser
 * boundary condition. Handles all four boundaries (top, bottom, left, right) for the local domain.
 */
void parallelSolver::initialize_sol() {
    solution = std::make_unique<matrix>(local_n + 2 - (first_rank || last_rank),
                                        param.grid_pts);

    // Bottom boundary condition
    if (first_rank) {
        param.boundary_point[1] = mesh_y[0];
        for (size_t i = 0; i < param.grid_pts; ++i) {
            param.boundary_point[0] = mesh_x[i];
            (*solution)(local_n, i) = param.boundary_condition.Eval();
        }
    }
    // Top boundary condition
    else if (last_rank) {
        param.boundary_point[1] = mesh_y[local_n - 1];
        for (size_t i = 0; i < param.grid_pts; ++i) {
            param.boundary_point[0] = mesh_x[i];
            (*solution)(0, i) = param.boundary_condition.Eval();
        }
    }

    // Left boundary condition
    param.boundary_point[0] = mesh_x[0];
    for (int i = 0; i < local_n; ++i) {
        param.boundary_point[1] = mesh_y[i];
        (*solution)(local_n - i, 0) = param.boundary_condition.Eval();
    }

    // Right boundary condition
    param.boundary_point[0] = mesh_x[param.grid_pts - 1];
    for (int i = 0; i < local_n; ++i) {
        param.boundary_point[1] = mesh_y[i];
        (*solution)(local_n - i, param.grid_pts - 1) =
            param.boundary_condition.Eval();
    }
}

/**
 * @brief Prints the current problem parameters from rank 0.
 * @details
 * Only the process with rank 0 prints the parameters using the overloaded operator<<.
 */
void parallelSolver::print_parameters() const {
    if (rank == 0) {
        std::cout << param << std::endl;
    }
}

/**
 * @brief Prints the mesh coordinates from rank 0.
 * @details
 * Gathers mesh_y data from all ranks and prints the full mesh in order, formatting each coordinate pair.
 */
void parallelSolver::print_mesh() const {
    if (rank > 0) {
        MPI_Send(mesh_y.data(), local_n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else {
        std::cout << "PRINTING THE MESH (from rank 0)" << std::endl;

        for (size_t i = size - 1; i > 0; --i) {
            int to_receive =
                (param.grid_pts / size) + (i < (param.grid_pts % size));

            std::vector<double> to_print(to_receive);

            MPI_Recv(to_print.data(), to_receive, MPI_DOUBLE, i, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int j = to_receive - 1; j >= 0; --j) {
                for (size_t k = 0; k < param.grid_pts; ++k) {
                    std::string s =
                        std::format("({:.2f}, {:.2f})", mesh_x[k], to_print[j]);
                    std::cout << std::setw(13) << s;
                }
                std::cout << std::endl;
            }
        }

        for (int j = local_n - 1; j >= 0; --j) {
            for (size_t k = 0; k < param.grid_pts; ++k) {
                std::string s =
                    std::format("({:.2f}, {:.2f})", mesh_x[k], mesh_y[j]);
                std::cout << std::setw(13) << s;
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;
    }
}

/**
 * @brief Prints the computed solution from rank 0.
 * @details
 * Gathers solution data from all ranks and prints the full solution grid in order.
 */
void parallelSolver::print_solution() const {
    if (rank > 0) {
        MPI_Send(solution->data() + !last_rank * param.grid_pts,
                 local_n * param.grid_pts, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else {
        std::cout << "PRINTING THE SOLUTION (from rank 0)" << std::endl;

        for (int i = size - 1; i > 0; --i) {
            int to_receive =
                (param.grid_pts / size) + (i < (param.grid_pts % size));

            std::vector<double> to_print(to_receive * param.grid_pts);

            MPI_Recv(to_print.data(), to_receive * param.grid_pts, MPI_DOUBLE,
                     i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (size_t j = 0; j < to_receive * param.grid_pts; ++j) {
                std::string s = std::format("{:.5e}", to_print[j]);
                std::cout << std::setw(13) << s;
                if ((j + 1) % param.grid_pts == 0) std::cout << std::endl;
            }
        }

        for (size_t i = 0; i < local_n; ++i) {
            for (size_t j = 0; j < param.grid_pts; ++j) {
                std::string s = std::format("{:.5e}", (*solution)(i + 1, j));
                std::cout << std::setw(13) << s;
            }

            std::cout << std::endl;
        }

        std::cout << std::endl;
    }
}

/**
 * @brief Reads parameters from a file and broadcasts them to all ranks.
 * @details
 * Rank 0 reads the JSON file and broadcasts all parameter values and expressions to other ranks.
 * If the file is not found, all ranks revert to default parameters.
 */
void parallelSolver::readParameters(std::string const& filename) {
    int len_f = 0, len_bc = 0, len_ex = 0;
    std::vector<char> f_buf, bc_buf, ex_buf;

    std::ifstream jfile;
    bool file_found = 1;
    if (rank == 0) {
        jfile.open(filename);
        if (!jfile) {
            std::cerr << "ERROR: parameters file " << filename
                      << " does not exist." << std::endl;
            std::cerr << "Reverting to default values." << std::endl;
            file_found = 0;
        }
    }

    MPI_Bcast(&file_found, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

    if (!file_found) {
        param.assign_defaults();
        return;
    }

    if (rank == 0) {
        nlohmann::json data = nlohmann::json::parse(jfile);

        param.max_it = data["max_iterations"].get<unsigned long long>();
        param.tolerance = data["tolerance"].get<double>();
        param.grid_pts =
            data["grid_points_in_each_direction"].get<unsigned long long>();
        param.exact_sol_given = data["exact_solution_given"].get<bool>();

        std::string fun_str = data["forcing_term"].get<std::string>();
        std::string bc_str = data["boundary_condition"].get<std::string>();
        std::string exact_str = data["exact_solution"].get<std::string>();
        jfile.close();

        len_f = fun_str.size() + 1;
        len_bc = bc_str.size() + 1;
        len_ex = exact_str.size() + 1;

        f_buf.resize(len_f);
        strncpy(f_buf.data(), fun_str.c_str(), len_f);
        f_buf[len_f - 1] = '\0';

        bc_buf.resize(len_bc);
        strncpy(bc_buf.data(), bc_str.c_str(), len_bc);
        bc_buf[len_bc - 1] = '\0';

        ex_buf.resize(len_ex);
        strncpy(ex_buf.data(), exact_str.c_str(), len_ex);
        ex_buf[len_ex - 1] = '\0';
    }

    MPI_Bcast(&len_f, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&len_bc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&len_ex, 1, MPI_INT, 0, MPI_COMM_WORLD);

    f_buf.resize(len_f);
    bc_buf.resize(len_bc);
    ex_buf.resize(len_ex);

    MPI_Bcast(f_buf.data(), len_f, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(bc_buf.data(), len_bc, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(ex_buf.data(), len_ex, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Bcast(&(param.max_it), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(param.tolerance), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(param.exact_sol_given), 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(param.grid_pts), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    mu::Parser forcing_term, boundary_condition, exact_solution;

    forcing_term.SetExpr(std::string{f_buf.begin(), f_buf.end() - 1});
    param.forcing_term = forcing_term;

    boundary_condition.SetExpr(std::string{bc_buf.begin(), bc_buf.end() - 1});
    param.boundary_condition = boundary_condition;

    exact_solution.SetExpr(std::string{ex_buf.begin(), ex_buf.end() - 1});
    param.exact_solution = exact_solution;

    param.define_variables();
}

/**
 * @brief Solves the Laplace equation using the Jacobi method in parallel.
 * @details
 * Allocates a new solution matrix and iteratively updates the solution using Jacobi sweeps.
 * Handles communication of boundary rows between neighboring MPI ranks and uses OpenMP for parallelization.
 * The convergence criterion is checked globally using MPI_Allreduce.
 */
void parallelSolver::solve() {
    std::unique_ptr<matrix> new_solution = std::make_unique<matrix>(
        local_n + 2 - (first_rank || last_rank), param.grid_pts);

#pragma omp parallel for
    // Left and right boundary conditions
    for (int i = 1; i <= local_n; ++i) {
        (*new_solution)(i, 0) = (*solution)(i, 0);
        (*new_solution)(i, param.grid_pts - 1) =
            (*solution)(i, param.grid_pts - 1);
    }

    matrix force;
    if (!(first_rank || last_rank) || local_n > 1) {
        force.resize(local_n - (first_rank || last_rank), param.grid_pts - 2);
    }

    // Building the force term
    for (size_t i = last_rank; i < local_n - first_rank; ++i) {
        for (size_t j = 1; j < param.grid_pts - 1; ++j) {
            param.force_point[0] = mesh_x[j];
            param.force_point[1] = mesh_y[local_n - 1 - i];
            force(i - last_rank, j - 1) =
                param.forcing_term.Eval() * mesh_size * mesh_size;
        }
    }

    bool not_finished = true;
    double conv_crit = 0.0;
    size_t it = 0;

    if (first_rank) {
#pragma omp parallel for
        // Bottom boundary condition
        for (size_t i = 0; i < param.grid_pts; ++i) {
            (*new_solution)(local_n, i) = (*solution)(local_n, i);
        }

        while (not_finished) {
            conv_crit = 0.0;

            MPI_Send(solution->data() + param.grid_pts, param.grid_pts,
                     MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(solution->data(), param.grid_pts, MPI_DOUBLE, 1, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#pragma omp parallel for reduction(+ : conv_crit)
            for (size_t i = 1; i < local_n; ++i) {
                for (size_t j = 1; j < param.grid_pts - 1; ++j) {
                    (*new_solution)(i, j) =
                        0.25 * ((*solution)(i - 1, j) + (*solution)(i + 1, j) +
                                (*solution)(i, j - 1) + (*solution)(i, j + 1) +
                                force(i - 1, j - 1));

                    conv_crit += ((*new_solution)(i, j) - (*solution)(i, j)) *
                                 ((*new_solution)(i, j) - (*solution)(i, j));
                }
            }

            std::swap(solution, new_solution);

            conv_crit = std::sqrt(mesh_size * conv_crit);
            ++it;

            not_finished = (it < param.max_it && conv_crit > param.tolerance);
            MPI_Allreduce(MPI_IN_PLACE, &not_finished, 1, MPI_CXX_BOOL, MPI_LOR,
                          MPI_COMM_WORLD);
        }
    }
    else if (last_rank) {
#pragma omp parallel for
        // Top boundary condition
        for (size_t i = 0; i < param.grid_pts; ++i) {
            (*new_solution)(0, i) = (*solution)(0, i);
        }

        if (rank % 2 == 0) {
            while (not_finished) {
                conv_crit = 0.0;

                MPI_Send(solution->data() + (local_n - 1) * param.grid_pts,
                         param.grid_pts, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD);
                MPI_Recv(solution->data() + local_n * param.grid_pts,
                         param.grid_pts, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#pragma omp parallel for reduction(+ : conv_crit)
                for (size_t i = 1; i < local_n; ++i) {
                    for (size_t j = 1; j < param.grid_pts - 1; ++j) {
                        (*new_solution)(i, j) =
                            0.25 *
                            ((*solution)(i - 1, j) + (*solution)(i + 1, j) +
                             (*solution)(i, j - 1) + (*solution)(i, j + 1) +
                             force(i - 1, j - 1));

                        conv_crit +=
                            ((*new_solution)(i, j) - (*solution)(i, j)) *
                            ((*new_solution)(i, j) - (*solution)(i, j));
                    }
                }

                std::swap(solution, new_solution);

                conv_crit = std::sqrt(mesh_size * conv_crit);
                ++it;

                not_finished =
                    (it < param.max_it && conv_crit > param.tolerance);
                MPI_Allreduce(MPI_IN_PLACE, &not_finished, 1, MPI_CXX_BOOL,
                              MPI_LOR, MPI_COMM_WORLD);
            }
        }
        else {
            while (not_finished) {
                conv_crit = 0.0;

                MPI_Recv(solution->data() + local_n * param.grid_pts,
                         param.grid_pts, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(solution->data() + (local_n - 1) * param.grid_pts,
                         param.grid_pts, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD);

#pragma omp parallel for reduction(+ : conv_crit)
                for (size_t i = 1; i < local_n; ++i) {
                    for (size_t j = 1; j < param.grid_pts - 1; ++j) {
                        (*new_solution)(i, j) =
                            0.25 *
                            ((*solution)(i - 1, j) + (*solution)(i + 1, j) +
                             (*solution)(i, j - 1) + (*solution)(i, j + 1) +
                             force(i - 1, j - 1));

                        conv_crit +=
                            ((*new_solution)(i, j) - (*solution)(i, j)) *
                            ((*new_solution)(i, j) - (*solution)(i, j));
                    }
                }

                std::swap(solution, new_solution);

                conv_crit = std::sqrt(mesh_size * conv_crit);
                ++it;

                not_finished =
                    (it < param.max_it && conv_crit > param.tolerance);
                MPI_Allreduce(MPI_IN_PLACE, &not_finished, 1, MPI_CXX_BOOL,
                              MPI_LOR, MPI_COMM_WORLD);
            }
        }
    }
    else {
        if (rank % 2 == 0) {
            while (not_finished) {
                conv_crit = 0.0;

                // Sending top row
                MPI_Send(solution->data() + param.grid_pts, param.grid_pts,
                         MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                // Receiving bottom row
                MPI_Recv(solution->data() + (local_n + 1) * param.grid_pts,
                         param.grid_pts, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // Sending bottom row
                MPI_Send(solution->data() + local_n * param.grid_pts,
                         param.grid_pts, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD);
                // Receiving top row
                MPI_Recv(solution->data(), param.grid_pts, MPI_DOUBLE, rank + 1,
                         0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#pragma omp parallel for reduction(+ : conv_crit)
                for (size_t i = 1; i <= local_n; ++i) {
                    for (size_t j = 1; j < param.grid_pts - 1; ++j) {
                        (*new_solution)(i, j) =
                            0.25 *
                            ((*solution)(i - 1, j) + (*solution)(i + 1, j) +
                             (*solution)(i, j - 1) + (*solution)(i, j + 1) +
                             force(i - 1, j - 1));

                        conv_crit +=
                            ((*new_solution)(i, j) - (*solution)(i, j)) *
                            ((*new_solution)(i, j) - (*solution)(i, j));
                    }
                }

                std::swap(solution, new_solution);

                conv_crit = std::sqrt(mesh_size * conv_crit);
                ++it;

                not_finished =
                    (it < param.max_it && conv_crit > param.tolerance);
                MPI_Allreduce(MPI_IN_PLACE, &not_finished, 1, MPI_CXX_BOOL,
                              MPI_LOR, MPI_COMM_WORLD);
            }
        }
        else {
            while (not_finished) {
                conv_crit = 0.0;

                // Receiving bottom row
                MPI_Recv(solution->data() + (local_n + 1) * param.grid_pts,
                         param.grid_pts, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // Sending top row
                MPI_Send(solution->data() + param.grid_pts, param.grid_pts,
                         MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                // Receiving top row
                MPI_Recv(solution->data(), param.grid_pts, MPI_DOUBLE, rank + 1,
                         0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // Sending bottom row
                MPI_Send(solution->data() + local_n * param.grid_pts,
                         param.grid_pts, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD);

#pragma omp parallel for reduction(+ : conv_crit)
                for (size_t i = 1; i <= local_n; ++i) {
                    for (size_t j = 1; j < param.grid_pts - 1; ++j) {
                        (*new_solution)(i, j) =
                            0.25 *
                            ((*solution)(i - 1, j) + (*solution)(i + 1, j) +
                             (*solution)(i, j - 1) + (*solution)(i, j + 1) +
                             force(i - 1, j - 1));

                        conv_crit +=
                            ((*new_solution)(i, j) - (*solution)(i, j)) *
                            ((*new_solution)(i, j) - (*solution)(i, j));
                    }
                }

                std::swap(solution, new_solution);

                conv_crit = std::sqrt(mesh_size * conv_crit);
                ++it;

                not_finished =
                    (it < param.max_it && conv_crit > param.tolerance);
                MPI_Allreduce(MPI_IN_PLACE, &not_finished, 1, MPI_CXX_BOOL,
                              MPI_LOR, MPI_COMM_WORLD);
            }
        }
    }

#ifdef DEBUG
    if (first_rank && it == param.max_it)
        std::cout
            << "Warning: parallel solver reached max number of iterations "
               "without "
               "satisfying the convergence criterion, stopping."
            << std::endl;
#endif
}

/**
 * @brief Sets the number of grid points per direction in parallel.
 * @param n Number of grid points.
 */
void parallelSolver::parallel_setpts(size_t n) { param.grid_pts = n; }

/**
 * @brief Exports the solution to a VTK file for visualization.
 * @details
 * Gathers the solution from all ranks and writes it to a VTK file in ASCII format.
 * Only rank 0 writes the file.
 * @param filename Output VTK file name.
 * @return True if export was successful, false otherwise.
 */
bool parallelSolver::export_vtk(std::string filename) const {
    if (rank > 0) {
        MPI_Send(solution->data() + !last_rank * param.grid_pts,
                 local_n * param.grid_pts, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else {
        std::ofstream vtkfile{filename};

        vtkfile << "# vtk DataFile Version 3.0" << std::endl;
        vtkfile << "Laplace Solver Solution" << std::endl;
        vtkfile << "ASCII" << std::endl;
        vtkfile << "DATASET STRUCTURED_POINTS" << std::endl;
        vtkfile << "DIMENSIONS " << param.grid_pts << " " << param.grid_pts
                << " "
                << "1" << std::endl;
        vtkfile << "ORIGIN 0.0 0.0 0.0" << std::endl;
        vtkfile << "SPACING " << mesh_size << " " << mesh_size << " " << "1.0"
                << std::endl;
        vtkfile << "POINT_DATA " << param.grid_pts * param.grid_pts
                << std::endl;
        vtkfile << "SCALARS solution_U double 1" << std::endl;
        vtkfile << "LOOKUP_TABLE default" << std::endl;

        for (int i = size - 1; i > 0; --i) {
            int to_receive =
                (param.grid_pts / size) + (i < (param.grid_pts % size));

            std::vector<double> to_print(to_receive * param.grid_pts);

            MPI_Recv(to_print.data(), to_receive * param.grid_pts, MPI_DOUBLE,
                     i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (size_t j = 0; j < to_receive * param.grid_pts; ++j) {
                std::string s = std::format("{:.5e}", to_print[j]);
                vtkfile << to_print[j] << std::endl;
            }
        }

        for (size_t i = 0; i < local_n; ++i) {
            for (size_t j = 0; j < param.grid_pts; ++j) {
                vtkfile << (*solution)(i + 1, j) << std::endl;
            }
        }
    }
    return true;
}

/**
 * @brief Computes and prints the error between the computed and exact solution.
 * @details
 * For each grid node, evaluates the exact solution and computes the squared difference
 * with the computed solution. The sum is reduced across all ranks and the square root is printed by rank 0.
 */
void parallelSolver::print_error() {
    double error = 0.0;

    for (size_t i = !last_rank + 1; i < local_n - !first_rank + 1; ++i) {
        for (size_t j = 1; j < param.grid_pts - 1; ++j) {
            param.exact_point[0] = mesh_x[j];
            param.exact_point[1] = mesh_y[i];
            double exact = param.exact_solution.Eval();
            error += (exact - (*solution)(j, i)) * (exact - (*solution)(j, i));
        }
    }

    MPI_Reduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "The square root of the sum of the squared differences "
                     "between the computeted solution and the provided exact "
                     "solution in the grid nodes is: "
                  << std::sqrt(error) << std::endl
                  << std::endl;
    }
}

}  // namespace laplace_solver
