[Table of Contents](README.md)

## Hadamard Product of Matrices

**Routine Name:** `hadamardProductOfMatrices`

**Author:** Kollin Murphy

**Language:** C

**Description/Purpose:** This routine will calculate the Hadamard product of two matrices. This function takes as input two matrices, `a` and `b`, and a third matrix `c` which is where the result will be stored. It also takes as input the number of rows and columns of the matrices.

It uses OpenMP to parallelize the computation.

#### Usage

The function accepts parameters as follows:

```c
hadamardProductOfMatrices(m, n, a, b, c)
```

Example usage:

```c
double a[4][4] = \{\{1, 2, 3, 4\}, \{5, 6, 7, 8\}, \{9, 10, 11, 12\}, \{13, 14, 15, 16\}\};
double b[4][4] = \{\{1, 2, 3, 4\}, \{5, 6, 7, 8\}, \{9, 10, 11, 12\}, \{13, 14, 15, 16\}\};
double c[4][4];
hadamardProductOfMatrices(4, 4, &a[0][0], &b[0][0], &c[0][0]);
for (int i = 0; i < 4; ++i) {
  for (int j = 0; j < 4; ++j) {
    printf("%f ", c[i][j]);
  }
  printf("\n");
}
```

Example output:

```
Time: 0.0000010000076145
1.000000 4.000000 9.000000 16.000000
25.000000 36.000000 49.000000 64.000000
81.000000 100.000000 121.000000 144.000000
169.000000 196.000000 225.000000 256.000000
```

#### Implementation

The function is implemented as follows:

```c
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
```
