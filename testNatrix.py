def rref(self, matrix):
    # Convert the matrix into a reduced row echelon form matrix
    lead = 0
    n = len(matrix)  # row_count
    m = len(matrix[0])  # column_count
    for r in range(n):
        if lead >= m:
            return matrix

        i = r
        while matrix[i][lead] == 0:
            i += 1
            if i == n:
                i = r
                lead += 1
                if m == lead:
                    return matrix

        matrix[i], matrix[r] = matrix[r], matrix[i]
        lead_value_r = matrix[r][lead]
        matrix[r] = [mrx / lead_value_r for mrx in matrix[r]]
        for i in range(n):
            if i != r:
                lead_value_i = matrix[i][lead]
                matrix[i] = [iv - lead_value_i * rv for rv, iv in zip(matrix[r], matrix[i])]
        lead += 1
    return matrix

def rref2(matrix):
    n, m = len(matrix), len(matrix[0])
    r, j = 0, 0
    while j < m and r < m:
        i = r
        while i < n:
            if matrix[i][j] != 0:
                matrix[i], matrix[r] = matrix[r], matrix[i]                         #החלפת שורות
                lead_value = matrix[r][j]
                matrix[r] = [matrix[r][k] / lead_value for k in range(m)]           #נרמול שורה מובילה
                for k in range(n):                                                  #איפוס השורות הבאות
                    lead_value_k = matrix[k][j]
                    matrix[k] = [kv - lead_value_k * rv for rv, kv in zip(matrix[r], matrix[k])] if k != r else matrix[r]
                r += 1
                break
            i += 1
        j += 1
    return matrix

def random_matrix(rows, cols):
    import random
    matrix = []
    for i in range(rows):
        matrix[i] = [random.randint(10) for j in range(cols)]
    return matrix

def inverse_matrix(matrix):
    n = len(matrix)
    for i in range(n):
        matrix[i] += [1 if i == j else 0 for j in range(n)]
    new_matrix = rref(matrix)
    for i in range(n):
        new_matrix[i] = new_matrix[i][n:]
    return new_matrix

def main_matrix():
    matrix1 = [[1,1,0,1,0,0],[1,0,1,0,1,0],[0,1,1,0,0,1]]
    matrix2 = [[1,1,0],[1,0,1],[0,1,1]]
    rref1_matrix = rref(matrix1)
    print(rref1_matrix)
    # print(matrix2)
    # rref2_matrix = rref2(matrix2)
    # print(rref2_matrix)
    inverse_matrix1 = inverse_matrix(matrix2)
    print(inverse_matrix1)
    # return print(rref1_matrix == rref2_matrix)

main_matrix()