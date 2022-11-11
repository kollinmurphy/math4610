from math import sqrt

def vectorAddition(v1, v2):
    """Returns the vector addition of v1 and v2"""
    if len(v1) != len(v2):
        raise ValueError("Vectors must be of same length")
    return [v1[i] + v2[i] for i in range(len(v1))]

def vectorSubtraction(v1, v2):
    """Returns the vector subtraction of v1 and v2"""
    if len(v1) != len(v2):
        raise ValueError("Vectors must be of same length")
    return [v1[i] - v2[i] for i in range(len(v1))]

def vectorScalarMultiplication(scalar, v):
    """Returns the scalar multiplication of scalar and v"""
    return [scalar * v[i] for i in range(len(v))]

def vectorL1Norm(v):
    """Returns the L1 norm of v"""
    return sum([abs(v[i]) for i in range(len(v))])

def vectorL2Norm(v):
    """Returns the L2 norm of v"""
    return sqrt(sum([v[i] * v[i] for i in range(len(v))]))

def vectorInfinityNorm(v):
    """Returns the infinity norm of v"""
    return max([abs(v[i]) for i in range(len(v))])

def vectorDotProduct(v1, v2):
    """Returns the dot product of v1 and v2"""
    if len(v1) != len(v2):
        raise ValueError("Vectors must be of same length")
    return sum([v1[i] * v2[i] for i in range(len(v1))])

def vectorCrossProduct(v1, v2):
    """Returns the cross product of v1 and v2"""
    if len(v1) != 3 or len(v2) != 3:
        raise ValueError("Vectors must be of length 3")
    return [v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]]

def vectorTripleProduct(v1, v2, v3):
    """Returns the triple product of v1, v2, and v3"""
    if len(v1) != 3 or len(v2) != 3 or len(v3) != 3:
        raise ValueError("Vectors must be of length 3")
    return vectorDotProduct(v1, vectorCrossProduct(v2, v3))

def actionOfMatrixOnVector(m, v):
    """Returns the action of matrix m on vector v"""
    if len(m[0]) != len(v):
        raise ValueError("Matrix and vector must be of appropriate length")
    return [vectorDotProduct(m[i], v) for i in range(len(m))]

def sumOfMatrices(m1, m2):
    """Returns the sum of matrices m1 and m2"""
    if len(m1) != len(m2):
        raise ValueError("Matrices must be of same length")
    return [vectorAddition(m1[i], m2[i]) for i in range(len(m1))]

def differenceOfMatrices(m1, m2):
    """Returns the difference of matrices m1 and m2"""
    if len(m1) != len(m2):
        raise ValueError("Matrices must be of same length")
    return [vectorSubtraction(m1[i], m2[i]) for i in range(len(m1))]

def col(matrix, i):
    """Returns the ith column of matrix"""
    if i >= len(matrix[0]) or i < 0:
        raise ValueError("Column index out of bounds")
    return [row[i] for row in matrix]

def productOfMatrices(m1, m2):
    """Returns the product of matrices m1 and m2"""
    if len(m1[0]) != len(m2):
        raise ValueError("Matrices must be of appropriate length")
    return [[vectorDotProduct(m1[i], col(m2, j)) for j in range(len(m2[0]))] for i in range(len(m1))]

def main():
    print(f"Vector Addition: vectorAddition([1, 2, 3], [4, 5, 6]) = {vectorAddition([1, 2, 3], [4, 5, 6])}")
    print(f"Vector Subtraction: vectorSubtraction([1, 2, 3], [4, 5, 6]) = {vectorSubtraction([1, 2, 3], [4, 5, 6])}")
    print(f"Vector Scalar Multiplication: vectorScalarMultiplication(2, [1, 2, 3]) = {vectorScalarMultiplication(2, [1, 2, 3])}")
    print(f"Vector Dot Product: vectorDotProduct([1, 2, 3], [4, 5, 6]) = {vectorDotProduct([1, 2, 3], [4, 5, 6])}")
    print(f"Sum of Matrices: sumOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2, 3], [4, 5, 6]]) = {sumOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2, 3], [4, 5, 6]])}")
    print(f"Difference of Matrices: differenceOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2, 3], [4, 5, 6]]) = {differenceOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2, 3], [4, 5, 6]])}")
    print(f"Product of Matrices: productOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2], [3, 4], [5, 6]]) = {productOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2], [3, 4], [5, 6]])}")
    print(f"Vector L1 Norm: vectorL1Norm([1, 2, 3]) = {vectorL1Norm([1, 2, 3])}")
    print(f"Vector L2 Norm: vectorL2Norm([1, 2, 3]) = {vectorL2Norm([1, 2, 3])}")
    print(f"Vector Infinity Norm: vectorInfinityNorm([1, 2, 3]) = {vectorInfinityNorm([1, 2, 3])}")
    print(f"Vector Cross Product: vectorCrossProduct([1, 2, 3], [4, 5, 6]) = {vectorCrossProduct([1, 2, 3], [4, 5, 6])}")
    print(f"Vector Triple Product: vectorTripleProduct([1, 2, 3], [4, 5, 6], [7, 8, 9]) = {vectorTripleProduct([1, 2, 3], [4, 5, 6], [7, 8, 9])}")
    print(f"Action of Matrix on Vector: actionOfMatrixOnVector([[1, 2, 3], [4, 5, 6]], [7, 8, 9]) = {actionOfMatrixOnVector([[1, 2, 3], [4, 5, 6]], [7, 8, 9])}")

main()
