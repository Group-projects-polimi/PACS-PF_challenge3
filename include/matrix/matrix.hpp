#ifndef Matrix_HPP
#define Matrix_HPP

#include <cassert>
#include <vector>

namespace algebra {

/**
 * @class matrix
 * @brief A simple matrix class for storing and manipulating 2D arrays of doubles.
 */
class matrix {
    using storage_type = std::vector<double>;

   private:
    size_t rows, columns; ///< Number of rows and columns in the matrix.
    storage_type internal_storage; ///< Internal storage for matrix elements.

   public:
    matrix() = default;

    /**
     * @brief Constructs a matrix with given number of rows and columns.
     * @param r Number of rows.
     * @param c Number of columns.
     */
    matrix(size_t r, size_t c) : rows{r}, columns{c} {
#ifdef DEBUG
        assert(r > 0 && c > 0 && "Error: wrong dimensions to matrix.");
#endif

        internal_storage.resize(r * c);
    }

    /**
     * @brief Const element access operator.
     * @param i Row index.
     * @param j Column index.
     * @return Value at (i, j).
     */
    double operator()(size_t i, size_t j) const;

    /**
     * @brief Element access operator.
     * @param i Row index.
     * @param j Column index.
     * @return Reference to value at (i, j).
     */
    double& operator()(size_t i, size_t j);

    /**
     * @brief Returns the number of rows.
     * @return Number of rows.
     */
    size_t get_rows() const;

    /**
     * @brief Returns the number of columns.
     * @return Number of columns.
     */
    size_t get_columns() const;

    /**
     * @brief Returns a pointer to the underlying data. Necessary for MPI.
     * @return Pointer to the data.
     */
    double* data();

    /**
     * @brief Returns the total number of elements.
     * @return Size of the matrix (rows * columns).
     */
    size_t size();

    /**
     * @brief Resizes the matrix to new dimensions.
     * @param r New number of rows.
     * @param c New number of columns.
     */
    void resize(size_t r, size_t c);
};

}  // namespace algebra

#endif
