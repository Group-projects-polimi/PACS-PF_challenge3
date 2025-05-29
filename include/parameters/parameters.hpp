#ifndef Parameters_HPP
#define Parameters_HPP

#include <array>

#include "muParser.h"

/**
 * @struct Parameters
 * @brief Stores configuration parameters for solving the Laplacian on a square domain using the Jacobi method.
 */
struct Parameters {
    /**
     * @brief Coordinates of a point where the forcing term is evaluated.
     */
    std::array<double, 2> force_point;
    /**
     * @brief Coordinates of a point on the boundary where the boundary condition is evaluated.
     */
    std::array<double, 2> boundary_point;
    /**
     * @brief Coordinates of a point where the exact solution is evaluated (if provided).
     */
    std::array<double, 2> exact_point;

    /**
     * @brief Maximum number of Jacobi iterations.
     */
    unsigned long long max_it{100};

    /**
     * @brief Convergence tolerance for the Jacobi method.
     */
    double tolerance{0.01};

    /**
     * @brief Number of grid points per spatial direction (the grid is grid_pts x grid_pts).
     */
    unsigned long long grid_pts{5};

    /**
     * @brief Whether the exact solution is provided.
     */
    bool exact_sol_given = true;

    /**
     * @brief muParser object for the forcing term f(x, y).
     */
    mu::Parser forcing_term;

    /**
     * @brief muParser object for the boundary condition g(x, y).
     */
    mu::Parser boundary_condition;

    /**
     * @brief muParser object for the exact solution u(x, y), if available.
     */
    mu::Parser exact_solution;

    /**
     * @brief Assigns default values to all parameters.
     */
    void assign_defaults();

    /**
     * @brief Defines variables for the muParser objects.
     */
    void define_variables();
};

/**
 * @brief Stream output operator for Parameters.
 * @param os Output stream.
 * @param params Parameters to print.
 * @return Reference to the output stream.
 */
std::ostream &operator<<(std::ostream &, const Parameters &);

#endif
