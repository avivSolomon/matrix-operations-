
"""
- מעבר בין בסיסים
- מציאת ערכים ווקטורים עצמיים
-מטריצה צמודה

לצמצם את הפונקציה של הדירוג והדטרמיננטה
להשלים את הטסטים לכל המתודות (בעקבות המעבר לצורת מחלקה- עצרתי בסכום מטריוצות)
"""


class Matrix:
    def __init__(self, matrix):
        """
        Constructor method for the Matrix class.

        Args:
            matrix (list): The input matrix.
        """
        self.matrix = matrix
        self.row = len(matrix)
        self.col = len(matrix[0])

    def rref(self, Matrix=None, Canonical_matrix=True):
        """
        Performs row reduction to its row echelon form (also known as reduced row echelon form) on the given matrix.

        Args:
            Matrix (list): The input matrix.
            Canonical_matrix (bool): Determines whether to convert the matrix to canonical form. Default is True.

        Returns:
            list: The row-reduced matrix.
        """
        matrix, n, m = (self.matrix, self.row, self.col) if Matrix is None or len(Matrix) == 0\
            else (Matrix, len(Matrix), len(Matrix[0]))
        mat = [row[:] for row in matrix]  # Create a copy of the matrix
        r, j = 0, 0
        while j < m and r < min(n, m):
            i = r
            while i < n:
                if mat[i][j] != 0:
                    mat[i], mat[r] = mat[r], mat[i]  # Swap rows
                    lead_value = mat[r][j]
                    mat[r][j:] = [mat[r][k] / lead_value for k in range(j, m)]  # Normalize leading row
                    for k in (range(n) if Canonical_matrix else range(r + 1, n)):  # Eliminate the column
                        lead_value_k = mat[k][j]
                        mat[k][j:] = [kv - lead_value_k * rv for rv, kv in
                                      zip(mat[r][j:], mat[k][j:])] if k != r else mat[r][j:]
                    r += 1
                    break
                i += 1
            j += 1
        return mat

    def inverse_matrix(self, Matrix=None):
        """
        Computes the inverse of the given matrix using row reduction.

        Args:
            matrix (list): The input matrix.

        Returns:
            list: The inverse matrix.
        """
        matrix = self.matrix if Matrix is None or len(Matrix) == 0 else Matrix
        n = self.row
        mat = [matrix[i][:] + [1 if i == j else 0 for j in range(n)] for i in range(n)]
        new_matrix = self.rref(mat)
        new_matrix = [new_matrix[i][n:] for i in range(n)]
        return new_matrix

    def transpose(self, Matrix=[]):
        """
        Computes the transpose of the given matrix.

        Args:
            matrix (list): The input matrix.

        Returns:
            list: The transposed matrix.
        """
        matrix = self.matrix if Matrix==[] else Matrix
        n,  m = self.row, self.col
        transpose_matrix = [[matrix[j][i] for j in range(m)] for i in range(n)]
        return transpose_matrix



    def matrix_sum(self, matrix2, Matrix=None):
        """
        Computes the sum of matrices A and B.

        Args:
            A (list): The first input matrix.
            B (list): The second input matrix.

        Returns:
            list: The sum matrix.
        """
        A, B = Matrix, matrix2
        n, m = self.row, self.col
        sum_matrix = [[]] * n
        for i in range(n):
            sum_matrix[i] = [A[i][j] + B[i][j] for j in range(m)]
        return sum_matrix

    def multiply_by_scaler(self, scaler, Matrix=None):
        """
        Multiplies the given matrix by a scalar.

        Args:
            scaler (float/int): The scalar value.
            matrix (list): The input matrix.

        Returns:
            list: The scaled matrix.
        """
        matrix = self.matrix if Matrix is None else Matrix
        return [[scaler * a for a in row] for row in matrix]

    def matrix_multiplication(self, matrix2, Matrix=None, mult_from_the_right=True):
        """
        Performs matrix multiplication between matrices A and B.

        Args:
            A (list): The first input matrix.
            B (list): The second input matrix.

        Returns:
            list: The resulting matrix.
        """
        matrix = self.matrix if Matrix is None else Matrix
        A, B = (matrix,matrix2) if mult_from_the_right else (matrix2, matrix)
        n, m, s = len(A), len(B[0]), len(A[0])
        multiplication = [[]] * n
        for i in range(n):
            multiplication[i] = [sum([A[i][k] * B[k][j] for k in range(s)]) for j in range(m)]
        return multiplication

    def determinant(self, Matrix=None):
        """
        Computes the determinant of the given matrix using row reduction.

        Args:
            matrix (list): The input matrix.

        Returns:
            float/int: The determinant value.
        """
        matrix = self.matrix if Matrix is None else Matrix
        n, r, j, det = self.row, 0, 0, 1
        mat = matrix.copy()
        while j < n and r < n:
            i = r
            while i < n:
                if matrix[i][j] != 0:
                    if i != r:
                        mat[i], mat[r], det = mat[r], mat[i], det * -1  # Swap rows
                    lead_value = matrix[r][j]
                    det *= matrix[r][j]
                    matrix[r] = [matrix[r][k] / lead_value for k in range(j, n)]  # Normalize leading row
                    for k in range(r + 1, n):  # Eliminate the column
                        lead_value_k = matrix[k][j]
                        matrix[k] = [kv - lead_value_k * rv for rv, kv in
                                     zip(matrix[r], matrix[k])] if k != r else matrix[r]
                    r += 1
                    break
                i += 1
            j += 1
        for i in range(n):
            det *= matrix[i][i]
        return det

    def get_characteristic_polynomial(self, Matrix=None):
        """
        Computes the coefficients of the characteristic polynomial.
        """
        matrix = self.matrix if Matrix is None else Matrix
        n = self.row
        coeffs = [1]
        for i in range(1, n + 1):
            sub_matrix = [matrix[j][:i] for j in range(i)]
            coeffs.append((-1) ** i * self.determinant(sub_matrix))
        return coeffs

    def eigenvalues(self, Matrix=None):
        """
        Computes the eigenvalues of the given square matrix.

        Args:
            matrix (list): The input square matrix.

        Returns:
            list: The eigenvalues of the matrix.
        """

        matrix = self.matrix if Matrix is None else Matrix
        matrix = self.matrix
        characteristic_coeffs = self.get_characteristic_polynomial(matrix)
        from numpy import roots as numpy_roots
        eigenvalues = numpy_roots(characteristic_coeffs)
        return eigenvalues.tolist()

    def trace(self, Matrix=None):
        """
        Computes the trace of the given matrix, which is the sum of the diagonal elements.

        Args:
            matrix (list): The input matrix.

        Returns:
            float/int: The trace value.
        """

        matrix = self.matrix if Matrix is None else Matrix
        return sum([matrix[i][i] for i in range(self.row)])

    def print_my_matrix(self):
        """
        Prints the matrix stored in the Matrix object.
        """
        print(self.matrix)


class Polynomial:
    def __init__(self, coeffs):
        self.coeffs = coeffs

    def evaluate(self, x):
        """
        Evaluates the polynomial at a given value of x.
        """
        result = 0
        for i in range(len(self.coeffs)):
            result += self.coeffs[i] * x**i
        return result

    def derivative(self):
        """
        Computes the derivative of the polynomial.
        """
        derivative_coeffs = [i * self.coeffs[i] for i in range(1, len(self.coeffs))]
        return Polynomial(derivative_coeffs)

    def newton_raphson(self, initial_guess, epsilon=1e-6, max_iterations=100):
        """
        Performs the Newton-Raphson method to find the roots of the polynomial.

        Args:
            initial_guess (float): The initial guess for the root.
            epsilon (float): The desired accuracy of the root. Defaults to 1e-6.
            max_iterations (int): The maximum number of iterations. Defaults to 100.

        Returns:
            float: The estimated root of the polynomial.
        """
        x = initial_guess
        iterations = 0

        while abs(self.evaluate(x)) > epsilon and iterations < max_iterations:
            x = x - self.evaluate(x) / self.derivative().evaluate(x)
            iterations += 1
        return x


"""----------------------------------------------------------------------"""

# Create a Matrix object
matrix_obj = Matrix([[1, 3], [4, 6]])

# # Test rref method
# mat = matrix_obj.rref(Canonical_matrix=True)
# print("Row-Reduced Matrix:")
# print(mat)

# # Test inverse_matrix method
# inverse_matrix = matrix_obj.inverse_matrix()
# print("Inverse Matrix:")
# print(inverse_matrix)

# # Test transpose method
# transpose_matrix = matrix_obj.transpose([[1, 2, 3], [4, 5, 6]])
# print("Transposed Matrix:")
# matrix_obj.print_my_matrix()

# Test matrix_sum method
sum_matrix = matrix_obj.matrix_sum([[1, 2], [3, 4]], [[5, 6], [7, 8]])
print("Sum Matrix:")
matrix_obj.print_my_matrix()

# # Test multiply_by_scaler method
# scaled_matrix = matrix_obj.multiply_by_scaler(2, [[1, 2], [3, 4]])
# print("Scaled Matrix:")
# matrix_obj.print_my_matrix()
#
# # Test matrix_multiplication method
# mult_matrix = matrix_obj.matrix_multiplication([[1, 2], [3, 4]], [[5, 6], [7, 8]])
# print("Multiplied Matrix:")
# matrix_obj.print_my_matrix()
#
# # Test determinant method
# det_value = matrix_obj.determinant([[1, 2], [3, 4]])
# print("Determinant:", det_value)
#
# # Test get_characteristic_polynomial method
# char_poly_coeffs = matrix_obj.get_characteristic_polynomial([[1, 2], [3, 4]])
# print("Characteristic Polynomial Coefficients:", char_poly_coeffs)
#
# # Test eigenvalues method
# eigenvalues = matrix_obj.eigenvalues([[1, 2], [3, 4]])
# print("Eigenvalues:", eigenvalues)
#
# # Test trace method
# trace_value = matrix_obj.trace([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
# print("Trace Value:", trace_value)
