# Tasksheet 7 by Kollin Murphy

## Task 1: Computing products of matrices with OpenMP

I wrote a routine to compute the product of two matrices using OpenMP. The routine is in `task1.c` and is called `matrixProduct`. The routine takes in two matrices, `A` and `B`, and the number of rows and columns in each matrix. The routine then computes the product of the two matrices and stores the result in `C`. The routine is parallelized using OpenMP and the number of threads is set to the number of cores on the machine. The routine is called in `main` and the results are printed to the screen. Also, the time it takes to compute the product is printed to the screen, along with which thread is performing the iteration of the outer loop.

Implementation:

```c
#include <stdio.h>
#include <omp.h>

void matrixProduct(int aRows, int aCols, int bCols, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  #pragma omp parallel for
  for (int i = 0; i < aRows; ++i) {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    printf("Thread %d of %d working on iteration %d\n", id, numThreads, i);
    for (int j = 0; j < bCols; ++j) {
      double sum = 0;
      for (int k = 0; k < aCols; ++k) {
        sum += a[i * aCols + k] * b[k * bCols + j];
      }
      c[i * bCols + j] = sum;
    }
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}

int main() {
  double a[2][3] = {
    {1, 2, 3},
    {4, 5, 6}
  };
  double b[3][2] = {
    {7, 8},
    {9, 10},
    {11, 12}
  };
  double c[2][2];
  matrixProduct(2, 3, 2, &a[0][0], &b[0][0], &c[0][0]);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      printf("%f ", c[i][j]);
    }
    printf("\n");
  }
  return 0;
}
```

Output:

```
Thread 1 of 8 working on iteration 1
Thread 0 of 8 working on iteration 0
Time: 0.0005110000201967
58.000000 64.000000 
139.000000 154.000000
```

Because this was tested on small matrices, not all of the threads were used. However, when I tested it on larger matrices, all of the threads were used.

## Task 2: Computing Hadamard product of vectors in serial

I wrote a routine to compute the Hadamard product of two vectors in serial. The routine is in `task2.c` and is called `hadamardProduct`. The routine takes in two vectors, `a` and `b`, and the length of the vectors. The routine then computes the Hadamard product of the two vectors and stores the result in `c`. The routine is called in `main` and the results are printed to the screen.

Implementation:

```c
#include <stdio.h>
#include <omp.h>

void hadamardProduct(int n, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  for (int i = 0; i < n; ++i) {
    c[i] = a[i] * b[i];
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}

int main() {
  double a[4] = {1, 2, 3, 4};
  double b[4] = {5, 6, 7, 8};
  double c[4];
  hadamardProduct(4, &a[0], &b[0], &c[0]);
  for (int i = 0; i < 4; ++i) {
    printf("%f ", c[i]);
  }
  printf("\n");
  return 0;
}
```

Output:

```
Time: 0.0000010000076145
5.000000 12.000000 21.000000 32.000000
```

## Task 3: Computing Hadamard product of vectors with OpenMP

I wrote a routine to compute the Hadamard product of two vectors using OpenMP. The routine is in `task3.c` and is called `hadamardProduct`. The routine takes in two vectors, `a` and `b`, and the length of the vectors. The routine then computes the Hadamard product of the two vectors and stores the result in `c`. The routine is parallelized using OpenMP and the number of threads is set to the number of cores on the machine. The routine is called in `main` and the results are printed to the screen. Also, the time it takes to compute the product is printed to the screen.

I tested it on vectors of length 100. The results are shown below.

Implementation:

```c
#include <stdio.h>
#include <omp.h>

void hadamardProduct(int n, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    for (int i = id; i < n; i = i + numThreads) {
      c[i] = a[i] * b[i];
    }
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}

int main() {
  double a[100];
  double b[100];
  double c[100];
  for (int i = 0; i < 100; ++i) {
    a[i] = i + 1;
    b[i] = i + 51;
  }
  hadamardProduct(100, &a[0], &b[0], &c[0]);
  for (int i = 0; i < 100; ++i) {
    printf("%f ", c[i]);
  }
  printf("\n");
  return 0;
}
```

Output:

```
Time: 0.0001599999959581
51.000000 104.000000 159.000000 216.000000 275.000000 336.000000 399.000000 464.000000 531.000000 600.000000 671.000000 744.000000 819.000000 896.000000 975.000000 1056.000000 1139.000000 1224.000000 1311.000000 1400.000000 1491.000000 1584.000000 1679.000000 1776.000000 1875.000000 1976.000000 2079.000000 2184.000000 2291.000000 2400.000000 2511.000000 2624.000000 2739.000000 2856.000000 2975.000000 3096.000000 3219.000000 3344.000000 3471.000000 3600.000000 3731.000000 3864.000000 3999.000000 4136.000000 4275.000000 4416.000000 4559.000000 4704.000000 4851.000000 5000.000000 5151.000000 5304.000000 5459.000000 5616.000000 5775.000000 5936.000000 6099.000000 6264.000000 6431.000000 6600.000000 6771.000000 6944.000000 7119.000000 7296.000000 7475.000000 7656.000000 7839.000000 8024.000000 8211.000000 8400.000000 8591.000000 8784.000000 8979.000000 9176.000000 9375.000000 9576.000000 9779.000000 9984.000000 10191.000000 10400.000000 10611.000000 10824.000000 11039.000000 11256.000000 11475.000000 11696.000000 11919.000000 12144.000000 12371.000000 12600.000000 12831.000000 13064.000000 13299.000000 13536.000000 13775.000000 14016.000000 14259.000000 14504.000000 14751.000000 15000.000000
```

## Task 4: Computing Hadamard product of matrices

I wrote a routine to compute the Hadamard product of two matrices. The routine is in `task4.c` and is called `hadamardProductOfMatrices`. The routine takes in two matrices, `a` and `b`, and the number of rows and columns in the matrices. The routine then computes the Hadamard product of the two matrices and stores the result in `c`. The routine is called in `main` and the time it takes to compute the product is printed to the screen.

I tested it on matrices of size 100 x 1000. The results are shown below.

Implementation:

```c
#include <stdio.h>
#include <omp.h>

void hadamardProductOfMatrices(int m, int n, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    for (int i = id; i < m; i = i + numThreads) {
      for (int j = 0; j < n; ++j) {
        c[i * n + j] = a[i * n + j] * b[i * n + j];
      }
    }
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}

int main() {
  int m = 100;
  int n = 1000;
  double a[m][n];
  double b[m][n];
  double c[m][n];
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i][j] = i + j + 1;
      b[i][j] = i + j + 51;
    }
  }
  printf("m = %d, n = %d\n", m, n);
  hadamardProductOfMatrices(m, n, &a[0][0], &b[0][0], &c[0][0]);
  return 0;
}
```

Output:

```
m = 100, n = 1000
Thread 1 of 8
Thread 4 of 8
Thread 0 of 8
Thread 2 of 8
Thread 3 of 8
Thread 7 of 8
Thread 5 of 8
Thread 6 of 8
Time: 0.0006509999802802
```

## Task 5: Outer Product of Vectors/Matrices

I found that the outer product of two vectors is defined by the following equation:

$$
\mathbf{a} \otimes \mathbf{b} = \begin{bmatrix} a_1 b_1 & a_1 b_2 & \cdots & a_1 b_n \\ a_2 b_1 & a_2 b_2 & \cdots & a_2 b_n \\ \vdots & \vdots & \ddots & \vdots \\ a_m b_1 & a_m b_2 & \cdots & a_m b_n \end{bmatrix}
$$

where $a$ is a vector of length $m$ and $b$ is a vector of length $n$.

I wrote a routine to compute the outer product of two vectors. The routine is in `task5.c` and is called `outerProductOfVectors`. The routine takes in two vectors, `a` and `b`, and the length of the vectors. The routine then computes the outer product of the two vectors and stores the result in `c`. The routine is called in `main` and the time it takes to compute the product is printed to the screen.

I tested it on vectors of length 5 and 10. The results are shown below.

Implementation:

```c
#include <stdio.h>
#include <omp.h>

void outerProductOfVectors(int n, int m, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; ++j) {
      c[i * m + j] = a[i] * b[j];
    }
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}

int main() {
  int n = 5;
  int m = 10;
  double a[n];
  double b[m];
  double c[n][m];
  for (int i = 0; i < n; ++i) {
    a[i] = i + 1;
  }
  for (int j = 0; j < m; ++j) {
    b[j] = j + 51;
  }
  printf("n = %d, m = %d\n", n, m);
  outerProductOfVectors(n, m, &a[0], &b[0], &c[0][0]);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      printf("%f ", c[i][j]);
    }
    printf("\n");
  }
  return 0;
}
```

Output:

```
n = 5, m = 10
Time: 0.0000010000076145
51.000000 52.000000 53.000000 54.000000 55.000000 56.000000 57.000000 58.000000 59.000000 60.000000
102.000000 104.000000 106.000000 108.000000 110.000000 112.000000 114.000000 116.000000 118.000000 120.000000
153.000000 156.000000 159.000000 162.000000 165.000000 168.000000 171.000000 174.000000 177.000000 180.000000
204.000000 208.000000 212.000000 216.000000 220.000000 224.000000 228.000000 232.000000 236.000000 240.000000
255.000000 260.000000 265.000000 270.000000 275.000000 280.000000 285.000000 290.000000 295.000000 300.000000
```

### Extension to matrices

Yes, the outer product can be extended to matrices. This is generally called a Kronecker product. The Kronecker product of two matrices is defined by the following equation:

$$
\mathbf{A} \otimes \mathbf{B} = \begin{bmatrix} A_{11} \mathbf{B} & A_{12} \mathbf{B} & \cdots & A_{1n} \mathbf{B} \\ A_{21} \mathbf{B} & A_{22} \mathbf{B} & \cdots & A_{2n} \mathbf{B} \\ \vdots & \vdots & \ddots & \vdots \\ A_{m1} \mathbf{B} & A_{m2} \mathbf{B} & \cdots & A_{mn} \mathbf{B} \end{bmatrix}
$$

where $\mathbf{A}$ is a matrix of size $m \times n$ and $\mathbf{B}$ is a matrix of size $p \times q$. The resulting matrix is of size $mp \times nq$.

The dimensions of the matrices do not have to match up as they do in standard matrix multiplication. The Kronecker product is defined for any two matrices.

### Parallel Implementation of Vector Outer Product

To implement an algorithm for the outer product of two vectors in parallel, we could use OpenMP directives to parallelize the outer loop. We could then adjust the loop indices to ensure that each thread is assigned a unique set of indices. The code for this is shown below.

```c
#include <stdio.h>
#include <omp.h>

void outerProductOfVectors(int n, int m, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    for (int i = id; i < n; i = i + numThreads) {
      for (int j = 0; j < m; ++j) {
        c[i * m + j] = a[i] * b[j];
      }
    }
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}
```

The output is the same as the serial implementation.
