"""
- מעבר בין בסיסים
- מציאת ערכים ווקטורים עצמיים
-מטריצה צמודה

צריך לתקן עוד את הפונקציה של מציאת המקדמים של הפולינום האופייני
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
        self.square_matrix = self.row == self.col
        self.canonical = self.rref(update_matrix=False)
        self.det = self.determinant(matrix) if self.square_matrix else None
        self.inverse = self.inverse_matrix(matrix, update_matrix=False) if self.det else None
        self.print = [[True, True, True, True], [0, None, 1]] # [matrix, details, canonical, inverse], [i, j, step]

    def set_print(self, matrix=True, details=False, canonical=False, inverse=False, i=0, j=None, step=1):
        """
        Prints the matrix stored in the Matrix object.

        Args:
            i (int): The starting index of rows. Default is 0.
            j (int): The ending index of rows. Default is None (all rows).
            step (int): The step size for rows. Default is 1.
        """
        self.print = [[matrix, details, canonical, inverse], [i, j, step]]

    def __str__(self):
        i, j, step = self.print[1]
        j = self.col if j is None else j
        mat = '\n'.join(['\t'.join([str(num) for num in row[i:j:step]])
                         for row in self.matrix])
        matrix = "\nThe Matrix:  \n" + mat + "\n"

        detail = "row: " + str(self.row) + "\n" + \
                 "col: " + str(self.col) + "\n" + \
                 "square matrix: " + str(self.square_matrix) + "\n"
        if self.square_matrix:
            detail += "The determinant: " + str(self.det) + "\n"

        con = '\n'.join(['\t'.join([str(num) for num in row])
                         for row in self.canonical])
        rref = "The canonical mode: \n" + con + "\n"

        if self.det:
            inverse = '\n'.join(['\t'.join([str(num) for num in row])
                                 for row in self.inverse])
            inv = "The inverse matrix: \n" + inverse + "\n"
        else:
            inv = ""

        s, bol = [matrix, detail, rref, inv], self.print[0]
        return '\n'.join([s[i] for i in range(4) if bol[i]])

    def rref(self, Matrix=None, Canonical_matrix=True, return_determinant=False, update_matrix=True):
        """
        Performs row reduction to its row echelon form (also known as reduced row echelon form) on the given matrix.

        Args:
            Matrix (list): The input matrix.
            Canonical_matrix (bool): Determines whether to convert the matrix to canonical form. Default is True.
            return_determinant (bool): Determines whether to return the determinant of the matrix. Default is False.
            update_matrix (bool): Determines whether to update the stored matrix with the row-reduced form. Default is True.

        Returns:
            list/float: The row-reduced matrix or the determinant value.
        """
        matrix, rows, columns = (Matrix, len(Matrix), len(Matrix[0])) if Matrix \
            else (self.matrix, self.row, self.col)
        mat = matrix if update_matrix else [row[:] for row in matrix]  # Create a copy of the matrix
        row_echelon, col, det = 0, 0, 1
        while col < columns and row_echelon < min(rows, columns):
            row = row_echelon
            while row < rows:
                if mat[row][col] != 0:
                    mat[row], mat[row_echelon], det = mat[row_echelon], mat[row], -det  # Swap rows
                    lead_value = mat[row_echelon][col]
                    det *= lead_value
                    mat[row_echelon][col:] = [mat[row_echelon][i] / lead_value for i in
                                              range(col, columns)]  # Normalize leading row
                    for k in (range(rows) if Canonical_matrix else range(row_echelon + 1, rows)):  # Eliminate the column
                        lead_value_k = mat[k][col]
                        mat[k][col:] = [kv - lead_value_k * rv for rv, kv in
                                        zip(mat[row_echelon][col:], mat[k][col:])] if k != row_echelon\
                            else mat[row_echelon][col:]
                    row_echelon += 1
                    break
                row += 1
            col += 1
        for i in range(rows):
            det *= mat[i][i]
        return det if return_determinant else mat

    def inverse_matrix(self, Matrix=None, update_matrix=True):
        """
        Computes the inverse of the given matrix using row reduction.

        Args:
            Matrix (list): The input matrix.
            update_matrix (bool): Determines whether to update the stored matrix with the inverse matrix. Default is True.

        Returns:
            list: The inverse matrix.
        """
        matrix, n = (Matrix, len(Matrix)) if Matrix else (self.matrix, self.row)
        mat = [matrix[i][:] + [1 if i == j else 0 for j in range(n)] for i in range(n)]
        new_matrix = self.rref(mat)
        new_matrix = [new_matrix[i][n:] for i in range(n)]
        if update_matrix:
            self.matrix = new_matrix
        return new_matrix

    def transpose(self, Matrix=None):
        """
        Computes the transpose of the given matrix.
        Args:
            Matrix (list): The input matrix.
        Returns:
            list: The transposed matrix.
        """
        matrix = Matrix if Matrix else self.matrix
        if Matrix:
            self.matrix = list(zip(*matrix))
            return self.matrix
        return list(zip(*matrix))

    def rotate_left(self, Matrix=None, update_matrix=True):
        matrix = Matrix if Matrix else self.matrix
        if update_matrix:
            self.matrix = list(zip(*matrix))[::-1]
            return self.matrix
        return list(zip(*matrix))[::-1]

    def matrix_sum(self, matrix2, Matrix=None, update_matrix=True):
        """
        Computes the sum of matrices A and B.

        Args:
            matrix2 (list): The second input matrix.
            Matrix (list): The first input matrix.
            update_matrix (bool): Determines whether to update the stored matrix with the sum matrix. Default is True.

        Returns:
            list: The sum matrix.
        """
        matrix, n, m = (self.matrix if update_matrix else self.matrix.copy(),self.row, self.col)\
            if Matrix is None else (Matrix, len(Matrix), len(Matrix[0]))
        if len(matrix2) != n and len(matrix2[0]) != m:
            raise ValueError("The matrices do not have the same dimension, so it is not possible to add them together")
        sum_matrix = [[matrix[i][j] + matrix2[i][j] for j in range(m)] for i in range(n)]
        return sum_matrix

    def multiply_by_scaler(self, scaler, Matrix=None, update_matrix=True):
        """
        Multiplies the given matrix by a scalar.

        Args:
            scaler (float/int): The scalar value.
            Matrix (list): The input matrix.
            update_matrix (bool): Determines whether to update the stored matrix with the scaled matrix. Default is True.

        Returns:
            list: The scaled matrix.
        """
        matrix = (self.matrix if update_matrix else self.matrix.copy()) if Matrix is None else Matrix
        matrix = [[scaler * a for a in row] for row in matrix]
        return matrix

    def matrix_multiplication(self, matrix2, Matrix=None, mult_from_the_right=True, update_matrix=True):
        """
        Performs matrix multiplication between matrices A and B.

        Args:
            matrix2 (list): The second input matrix.
            Matrix (list): The first input matrix.
            mult_from_the_right (bool): Determines whether to multiply A from the right side of B. Default is True.
            update_matrix (bool): Determines whether to update the stored matrix with the resulting matrix. Default is True.

        Returns:
            list: The resulting matrix.
        """
        if Matrix is None:
            A, B = self.matrix if update_matrix else self.matrix.copy(), matrix2
        else:
            A, B = matrix2, Matrix
        if not mult_from_the_right: A, B = B, A

        if len(A[0]) != len(B):
            raise ValueError("The number of columns of the first matrix does not match the number of rows of the second,"
                             " so they cannot be multiplied")

        n, m, s = len(A), len(B[0]), len(A[0])
        multiplication = [[]] * n
        for i in range(n):
            multiplication[i] = [sum([A[i][k] * B[k][j] for k in range(s)]) for j in range(m)]
        return multiplication


    def determinant(self, Matrix=None):
        """
        Computes the determinant of the given matrix using row reduction.

        Args:
            Matrix (list): The input matrix.

        Returns:
            float/int: The determinant value.
        """
        matrix = Matrix if Matrix else self.matrix
        if len(matrix) != len(matrix[0]):
            raise ValueError("The matrix is not square, so the determinant cannot be computed for it")
        return self.rref(Matrix=matrix, return_determinant=True, update_matrix=False)

    def recursive_determinant(self,det=0, Matrix=None):
        """
        Computes the determinant of the given matrix recursively.

        Args:
            det (float/int): The current determinant value. Default is 0.
            Matrix (list): The input matrix.

        Returns:
            float/int: The determinant value.
        """
        n = len(Matrix)
        if n <= 2:
            det += Matrix[1][1]*Matrix[2][2]-Matrix[1][2]*Matrix[2][1]
            return det
        else:
            for i in range(n):
                open_row = Matrix.pop(i)
                for j in range(len(open_row)):
                    new_matrix = self.transpose(Matrix, update_matrix=False)
                    new_matrix.pop(j)
                    new_matrix = self.transpose(new_matrix)
                    det += (-1)**(i+j)*(open_row[j])*self.recursive_determinant(det, new_matrix)
                    return det
    def get_characteristic_polynomial(self, Matrix=None):
        """
        Computes the coefficients of the characteristic polynomial.

        Args:
            Matrix (list): The input matrix.

        Returns:
            list: The coefficients of the characteristic polynomial.
        """
        if Matrix is None:
            matrix, n = self.matrix, self.row
        else:
            matrix, n = Matrix, len(Matrix)

        coeffs = [1, -self.trace(matrix)]
        for i in range(2, n + 1):
            sub_matrix = [matrix[j][:i] for j in range(i)]
            det = self.determinant(sub_matrix)
            coeffs.append((-1) ** (n - i) * det)
        return coeffs

    def eigenvalues(self, Matrix=None):
        """
        Computes the eigenvalues of the given square matrix.

        Args:
            Matrix (list): The input square matrix.

        Returns:
            list: The eigenvalues of the matrix.
        """
        # matrix = self.matrix if Matrix is None else Matrix
        # characteristic_coeffs = self.get_characteristic_polynomial(matrix)
        # from numpy import roots as numpy_roots
        # eigenvalues = numpy_roots(characteristic_coeffs)
        # return eigenvalues.tolist()

    def reshape(self, r :int = 1, c :int = None, Matrix :list[list] = None) -> list[list[int]]:
        matrix, n, m = (Matrix, len(Matrix), len(Matrix[0])) if Matrix else (self.matrix, self.row, self.col)
        flag = True if c is None and r else False
        if r * c != m * n:
            return matrix
        arr = []
        for row in matrix:
            arr += row
        return arr if flag else [arr[c * i: c * (i + 1)] for i in range(r)]

    def trace(self, Matrix=None, return_list=False):
        """
        Computes the trace of the given matrix, which is the sum of the diagonal elements.

        Args:
            Matrix (list): The input matrix.
            return_list (bool): Determines whether to return a list of diagonal elements. Default is False.

        Returns:
            float/int/list: The trace value or a list of diagonal elements.
        """
        matrix, n = (Matrix, len(Matrix)) if Matrix else self.matrix, self.row
        trace_list = [matrix[i][i] for i in range(n)]
        return trace_list if return_list else sum(trace_list)

    def subset(self, lst, allSet=[], i=0, cur=[]):
        """
        Finds all possible subsets of a given list.

        Args:
            lst (list): The input list.
            allSet (list): The list to store all subsets. Default is an empty list.
            i (int): The current index. Default is 0.
            cur (list): The current subset. Default is an empty list.

        Returns:
            list: The list of all subsets.
        """
        if i >= len(lst):
            allSet.append(cur.copy())
            return allSet

        cur.append(lst[i])
        allSet = self.subset(lst, allSet, i + 1, cur)
        cur.pop()
        allSet = self.subset(lst, allSet, i + 1, cur)
        return allSet

