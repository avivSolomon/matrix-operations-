
""""
בניית מחלקה למטריצות
פעולות שיהיה ניתן לבצע:
- דירוג וקינון - היפוך
- שיחלוף
- כפל מטריצה בסקלר
- חיבור מטריצות
- כפל מטריצות

- מעבר בין בסיסים
- מציאת ערכים ווקטורים עצמיים
- דטרמיננטה
- עקבה

-מטריצה צמודה

"""


class Matrix:
    def __init__(self,matrix):
        self.matrix = matrix
        self.row = len(matrix)
        self.col = len(matrix[0])

    def rref(self, matrix, Canonical_matrix=True):
        n, m = len(matrix), len(matrix[0])
        r, j = 0, 0
        while j < m and r < min(n, m):
            i = r
            while i < n:
                if matrix[i][j] != 0:
                    matrix[i], matrix[r] = matrix[r], matrix[i]                         #החלפת שורות
                    lead_value = matrix[r][j]
                    matrix[r] = [matrix[r][k] / lead_value for k in range(m)]           #נרמול שורה מובילה
                    for k in range(0 if Canonical_matrix else r+1, n):                                    # איפוס העמודה
                        lead_value_k = matrix[k][j]
                        matrix[k] = [kv - lead_value_k * rv for rv, kv in zip(matrix[r], matrix[k])] if k != r else matrix[r]
                    r += 1
                    break
                i += 1
            j += 1
        return matrix


    def inverse_matrix(self, matrix):
        n = len(matrix)
        for i in range(n):
            matrix[i] += [1 if i == j else 0 for j in range(n)]
        new_matrix = rref(matrix)
        for i in range(n):
            new_matrix[i] = new_matrix[i][n:]
        return new_matrix

    def transpose(self, matrix):
        m, n = len(matrix), len(matrix[0])
        transpose_matrix = [[]]*n
        for i in range(n):
            transpose_matrix[i] = ([matrix[j][i] for j in range(m)])
        return transpose_matrix

    def Matrix_sum(self, A, B):
        n, m = len(A), len(A[0])
        sum_matrix = [[]] * n
        for i in range(n):
            sum_matrix[i] = [A[i][j] + B[i][j] for j in range(m)]
        return sum_matrix

    def multiplyBy_scaler(self, scaler, matrix):
        n = len(A)
        new_matrix = [[]]*n
        for i in range(n):
            new_matrix[i] = [scaler * a for a in matrix[i]]
        return new_matrix

    def Matrix_multiplication(self, A, B):
        n, m, s = len(A), len(B[0]), len(A[0])
        multiplication = [[]] * n
        for i in range(n):
            multiplication[i] = [sum([A[i][k] * B[k][j] for k in range(s)]) for j in range(m)]
        return multiplication

    def print_my_matrix(self):
            return print(self.matrix)



# if __name__ == '__main__':
#     # Example usage
#     matrix = [[4, 1, 1, 0], [6, 4, 0, 1]]
#     # rref(matrix)
#     # myMatrix = random_matrix()
#     myMatrix = matrix
#     print("~~~~~~~~~~~before~~~~~~~~~~~~~~~~~")
#     for r0w in range(len(myMatrix)):
#         print(myMatrix[rew])
#     rref(myMatrix)
#     print("~~~~~~~~~~~after~~~~~~~~~~~~~~~~~")
#     for row in range(len(myMatrix)):
#         print(myMatrix[rew])

