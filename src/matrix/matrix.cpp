#include "matrix.hpp"

namespace algebra {

/**
 * @brief Const element access operator.
 * @details Returns the value at position (i, j) in the matrix.
 * Performs bounds checking in debug mode.
 * @param i Row index (zero-based).
 * @param j Column index (zero-based).
 * @return Value at (i, j).
 */
double matrix::operator()(size_t i, size_t j) const {
#ifdef DEBUG
    assert(i >= 0 && i < rows && j >= 0 && j < columns &&
           "Error: wrong indexes to matrix.");
#endif

    return internal_storage[i * columns + j];
}

/**
 * @brief Element access operator.
 * @details Returns a reference to the value at position (i, j) in the matrix.
 * Performs bounds checking in debug mode.
 * @param i Row index (zero-based).
 * @param j Column index (zero-based).
 * @return Reference to value at (i, j).
 */
double& matrix::operator()(size_t i, size_t j) {
#ifdef DEBUG
    assert(i >= 0 && i < rows && j >= 0 && j < columns &&
           "Error: wrong indexes to matrix.");
#endif

    return internal_storage[i * columns + j];
}

/**
 * @brief Returns a pointer to the underlying data array.
 * @details Useful for interoperability with libraries (e.g., MPI).
 * @return Pointer to the first element of the internal storage.
 */
double* matrix::data() { return internal_storage.data(); }

/**
 * @brief Returns the total number of elements in the matrix.
 * @details Equivalent to rows * columns.
 * @return Number of elements.
 */
size_t matrix::size() { return internal_storage.size(); };

/**
 * @brief Resizes the matrix to new dimensions.
 * @details Updates the number of rows and columns, and resizes the internal storage.
 * Performs bounds checking in debug mode.
 * @param r New number of rows.
 * @param c New number of columns.
 */
void matrix::resize(size_t r, size_t c) {
#ifdef DEBUG
    assert(r > 0 && c > 0 && "Error: wrong dimensions to resize method.");
#endif

    rows = r;
    columns = c;
    internal_storage.resize(r * c);
}

/**
 * @brief Returns the number of rows in the matrix.
 * @return Number of rows.
 */
size_t matrix::get_rows() const { return rows; }

/**
 * @brief Returns the number of columns in the matrix.
 * @return Number of columns.
 */
size_t matrix::get_columns() const { return columns; }

}  // namespace algebra
