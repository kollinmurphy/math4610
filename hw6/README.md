# Tasksheet 6 by Kollin Murphy

## Task 1: Computing an Approximation of Pi

To calculate an approximation of $\pi$, I decided to use the following integral:
$$
\int \frac{1}{\sqrt{1-x^2}}dx = \arcsin(x)+C
$$
 By taking the definite integral from $0$ to $1$, we are able to get a multiple of $\pi$.
$$
\begin{align*}
\int_0^{\frac{1}{2}} \frac{1}{\sqrt{1-x^2}}dx &= \arcsin(\frac{1}{2})-\arcsin(0) \\
&= \frac{\pi}{6}-0 \\
&= \frac{\pi}{6}
\end{align*}
$$
Therefore, we can use the following to generate our approximation:
$$
6\int_0^{\frac{1}{2}} \frac{1}{\sqrt{1-x^2}}dx=\pi
$$

## Task 2: Coding an Approximation of Pi

To code an approximation of $\pi$, I used the trapezoid rule to approximate the integral. I used the following code to do this:

```c
#include <stdio.h>
#include <math.h>

double f(double x) {
  return 1 / sqrtf(1 - (x * x));
}

double trapezoidRule(double (*f)(double), double a, double b, int n) {
  double h = (b - a) / n;
  double sum = 0.5 * (f(a) + f(b));
  for (int i = 1; i < n; ++i) {
    sum += f(a + i * h);
  }
  return sum * h;
}

int main() {
  double pi = 6 * trapezoidRule(f, 0, 0.5, 100000000);
  printf("pi = %.16f\n", pi);
  return 0;
}
```

This code is serial. It uses a total of $100000000$ iterations to approximate $\pi$. The output of this code is:

```
pi = 3.1415926573648885
```

It ran in $0.201$ seconds on my machine as timed by using the `time` command.

## Task 3: Parallelizing an Approximation of Pi

I added OpenMP directives to parallelize the code. This is my implementation:

```c
#include <stdio.h>
#include <math.h>
#include <omp.h>

double f(double x) {
  return 1 / sqrtf(1 - (x * x));
}

static long num_steps = 100000000;
#define NUM_THREADS 4

double trapezoidRule(double (*f)(double), double a, double b, int n) {
  double h = (b - a) / n;
  double sum[NUM_THREADS];
  int actualNumThreads = 0;
  double accumulator = 0.5 * (f(a) + f(b));

  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    if (id == 0) actualNumThreads = numThreads;

    for (int i = 1 + id; i < n; i = i + numThreads) {
      sum[id] += f(a + i * h);
    }
  }

  for (int i = 0; i < actualNumThreads; ++i) {
    accumulator += sum[i];
  }
  return accumulator * h;
}

int main() {
  double pi = 6 * trapezoidRule(f, 0, 0.5, num_steps);
  printf("pi = %.16f\n", pi);
  return 0;
}
```

As before, I got the same output `pi = 3.1415926573648885`. With the compiler optimizations and set to run on 4 threads, it ran in $0.174$ seconds, a little bit faster than the serial code but definitely not as much of a decrease as I expected. I assume that this is due to increased overhead by setting up and taking down the threads. The calculations done here just don't take long enough to start seeing massive benefits from parallelization.

## Task 4: Approximating the value of $e$



## Task 5: Linear Algebra Algorithms

I implemented each algorithm. Below is the implementation of each function.

```python
from math import sqrt

def sumOfVectors(v1, v2):
    """Returns the sum of vector v1 and vector v2"""
    if len(v1) != len(v2):
        raise ValueError("Vectors must be of same length")
    return [v1[i] + v2[i] for i in range(len(v1))]

def differenceOfVectors(v1, v2):
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

def matrixVectorProduct(m, v):
    """Returns the action of matrix m on vector v"""
    if len(m[0]) != len(v):
        raise ValueError("Matrix and vector must be of appropriate length")
    return [vectorDotProduct(m[i], v) for i in range(len(m))]

def sumOfMatrices(m1, m2):
    """Returns the sum of matrices m1 and m2"""
    if len(m1) != len(m2):
        raise ValueError("Matrices must be of same length")
    return [sumOfVectors(m1[i], m2[i]) for i in range(len(m1))]

def differenceOfMatrices(m1, m2):
    """Returns the difference of matrices m1 and m2"""
    if len(m1) != len(m2):
        raise ValueError("Matrices must be of same length")
    return [differenceOfVectors(m1[i], m2[i]) for i in range(len(m1))]

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
```

The following is the test function I used to test the algorithms.

```python
print(f"Vector Addition: sumOfVectors([1, 2, 3], [4, 5, 6]) = {sumOfVectors([1, 2, 3], [4, 5, 6])}")
print(f"Vector Subtraction: differenceOfVectors([1, 2, 3], [4, 5, 6]) = {differenceOfVectors([1, 2, 3], [4, 5, 6])}")
print(f"Vector Scalar Multiplication: vectorScalarMultiplication(2, [1, 2, 3]) = {vectorScalarMultiplication(2, [1, 2, 3])}")
print(f"Vector Dot Product: vectorDotProduct([1, 2, 3], [4, 5, 6]) = {vectorDotProduct([1, 2, 3], [4, 5, 6])}")
print(f"Sum of Matrices: sumOfMatrices([[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]) = {sumOfMatrices([[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]])}")
print(f"Difference of Matrices: differenceOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2, 3], [4, 5, 6]]) = {differenceOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2, 3], [4, 5, 6]])}")
print(f"Product of Matrices: productOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2], [3, 4], [5, 6]]) = {productOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2], [3, 4], [5, 6]])}")
print(f"Vector L1 Norm: vectorL1Norm([1, 2, 3]) = {vectorL1Norm([1, 2, 3])}")
print(f"Vector L2 Norm: vectorL2Norm([1, 2, 3]) = {vectorL2Norm([1, 2, 3])}")
print(f"Vector Infinity Norm: vectorInfinityNorm([1, 2, 3]) = {vectorInfinityNorm([1, 2, 3])}")
print(f"Vector Cross Product: vectorCrossProduct([1, 2, 3], [4, 5, 6]) = {vectorCrossProduct([1, 2, 3], [4, 5, 6])}")
print(f"Vector Triple Product: vectorTripleProduct([1, 2, 3], [4, 5, 6], [7, 8, 9]) = {vectorTripleProduct([1, 2, 3], [4, 5, 6], [7, 8, 9])}")
print(f"Action of Matrix on Vector: matrixVectorProduct([[1, 2, 3], [4, 5, 6]], [7, 8, 9]) = {matrixVectorProduct([[1, 2, 3], [4, 5, 6]], [7, 8, 9])}")
```

The following is the output of the test function.

```
Vector Addition: sumOfVectors([1, 2, 3], [4, 5, 6]) = [5, 7, 9]
Vector Subtraction: differenceOfVectors([1, 2, 3], [4, 5, 6]) = [-3, -3, -3]
Vector Scalar Multiplication: vectorScalarMultiplication(2, [1, 2, 3]) = [2, 4, 6]
Vector Dot Product: vectorDotProduct([1, 2, 3], [4, 5, 6]) = 32
Sum of Matrices: sumOfMatrices([[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]) = [[8, 10, 12], [14, 16, 18]]
Difference of Matrices: differenceOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2, 3], [4, 5, 6]]) = [[0, 0, 0], [0, 0, 0]]
Product of Matrices: productOfMatrices([[1, 2, 3], [4, 5, 6]], [[1, 2], [3, 4], [5, 6]]) = [[22, 28], [49, 64]]
Vector L1 Norm: vectorL1Norm([1, 2, 3]) = 6
Vector L2 Norm: vectorL2Norm([1, 2, 3]) = 3.7416573867739413
Vector Infinity Norm: vectorInfinityNorm([1, 2, 3]) = 3
Vector Cross Product: vectorCrossProduct([1, 2, 3], [4, 5, 6]) = [-3, 6, -3]
Vector Triple Product: vectorTripleProduct([1, 2, 3], [4, 5, 6], [7, 8, 9]) = 0
Action of Matrix on Vector: matrixVectorProduct([[1, 2, 3], [4, 5, 6]], [7, 8, 9]) = [50, 122]
```

You can find further documentation for each of these functions at https://kollinmurphy.github.io/math4610/docs/.
