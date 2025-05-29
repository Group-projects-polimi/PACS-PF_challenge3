#include "parameters.hpp"

/**
 * @brief Defines variables for the muParser objects.
 * @details
 * For each spatial dimension, defines a variable (x1, x2, ...) in the muParser objects
 * for the forcing term, boundary condition, and exact solution. The variables are linked
 * to the corresponding entries in force_point, boundary_point, and exact_point.
 */
void Parameters::define_variables() {
    std::string var{"x0"};
    for (size_t i = 0; i < force_point.size(); ++i) {
        var.erase(var.length() - (i + 1) / 10 - 1);
        var += std::to_string(i + 1);

        forcing_term.DefineVar(var, &force_point[i]);
        boundary_condition.DefineVar(var, &boundary_point[i]);
        exact_solution.DefineVar(var, &exact_point[i]);
    }
}

/**
 * @brief Assigns default values to all parameters.
 * @details
 * Sets default values for the maximum number of iterations, tolerance, grid points,
 * and whether the exact solution is given. Also sets default expressions for the
 * forcing term, boundary condition, and exact solution, and defines the variables
 * in the muParser objects.
 */
void Parameters::assign_defaults() {
    max_it = 100;
    tolerance = 0.01;
    grid_pts = 5;
    exact_sol_given = true;

    forcing_term.SetExpr("8*_pi^2*sin(2*_pi*x1)*sin(2*_pi*x2)");
    boundary_condition.SetExpr("0");
    exact_solution.SetExpr("sin(2*_pi*x1)*sin(2*_pi*x2)");

    define_variables();
}

/**
 * @brief Stream output operator for Parameters.
 * @details
 * Prints the values of all parameters, including the expressions for the forcing term,
 * boundary condition, and exact solution, to the provided output stream.
 * @param out Output stream.
 * @param p Parameters object to print.
 * @return Reference to the output stream.
 */
std::ostream &operator<<(std::ostream &out, Parameters const &p) {
    out << "PARAMETERS' VALUES:" << std::endl;
    out << "Max iterations = " << p.max_it << std::endl;
    out << "Tolerance = " << p.tolerance << std::endl;
    out << "Grid points = " << p.grid_pts << std::endl;
    out << "Exact solution given = " << p.exact_sol_given << std::endl;

    out << "Forcing term = " << p.forcing_term.GetExpr() << std::endl;
    out << "Boundary condition = " << p.boundary_condition.GetExpr()
        << std::endl;
    out << "Exact solution = " << p.exact_solution.GetExpr() << std::endl;

    return out;
}
