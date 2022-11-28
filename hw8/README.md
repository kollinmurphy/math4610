# Tasksheet 9 by Kollin Murphy

## Task 1: Kronecker Product of Two Matrices

I implemented the Kronecker product of two matrices in the function `computeKroneckerProduct` in `task1.c`. The function takes in three arrays of doubles, two for the two matrices and one for the result. The function also takes in the number of rows and columns for the input matrices. The function then computes the Kronecker product of the two matrices and stores the result in the result array. The function returns void because the result is stored in the result array.

I implemented it using OpenMP directives to parallelize the computation. It also prints out the threads used and the time it took to compute the Kronecker product.

Implementation:

```c
#include <stdio.h>
#include <omp.h>

void computeKroneckerProduct(double *A, double *B, double *C, int aRows, int aCols, int bRows, int bCols)
{
  int i, j, k, l;
  double start = omp_get_wtime();
#pragma omp parallel for
  for (i = 0; i < aRows; i++)
  {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    printf("Thread %d of %d\n", id, numThreads);
    for (j = 0; j < aCols; j++)
    {
      for (k = 0; k < bRows; k++)
      {
        for (l = 0; l < bCols; l++)
        {
          C[(i * bRows + k) * (aCols * bCols) + (j * bCols + l)] = A[i * aCols + j] * B[k * bCols + l];
        }
      }
    }
  }
  double end = omp_get_wtime();
  printf("Time: %f\n", end - start);
}
```

I tested the function with the following code:

```c
int main() {
  int aRows = 2;
  int aCols = 2;
  int bRows = 2;
  int bCols = 2;
  double A[aRows * aCols];
  double B[bRows * bCols];
  double C[aRows * bRows * aCols * bCols];
  A[0] = 1;
  A[1] = 2;
  A[2] = 3;
  A[3] = 4;
  B[0] = 0;
  B[1] = 5;
  B[2] = 6;
  B[3] = 7;
  computeKroneckerProduct(A, B, C, aRows, aCols, bRows, bCols);

  for (int i = 0; i < aRows * bRows; i++)
  {
    for (int j = 0; j < aCols * bCols; j++)
    {
      printf("%8.2f ", C[i * aCols * bCols + j]);
    }
    printf("\n");
  }
  return 0;
}
```

Output:

```bash
Thread 1 of 8
Thread 0 of 8
Time: 0.000683
    0.00     5.00     0.00    10.00 
    6.00     7.00    12.00    14.00 
    0.00    15.00     0.00    20.00 
   18.00    21.00    24.00    28.00
```

The function was able to successfully compute the Kronecker product of the two matrices. The output is correct.

## Task 2: Serial Power Method

I implemented the power method in serial in C. The function is called `computePowerMethod` and is located in `task2.c`. The function takes in a square matrix, a vector, and the number of rows and columns in the matrix. It also takes in a maximum number of iterations and a tolerance. The function then uses the power method on the matrix to approximate its largest eigenvalue and eigenvector. The function returns the approximate eigenvalue and stores the approximate eigenvector in the vector array.

To simplify the code I extracted some of the code into helper functions. These are:

```c
void matrixVectorProduct(int n, int m, double *matrix, double *vector, double *result)
{
  double s;
  for (int i = 0; i < n; i++)
  {
    result[i] = 0;
    for (int j = 0; j < m; j++)
    {
      result[i] += matrix[i * n + j] * vector[j];
    }
  }
}

double l2Norm(int n, double vector[n])
{
  double sum = 0;
  for (int i = 0; i < n; i++)
  {
    sum += vector[i] * vector[i];
  }
  return sqrtf(sum);
}

double dotProduct(int n, double vector1[n], double vector2[n])
{
  double sum = 0;
  for (int i = 0; i < n; i++)
  {
    sum += vector1[i] * vector2[i];
  }
  return sum;
}
```

The implementation of the power method is as follows:

```c
double computePowerMethod(double *A, int n, double *v0, int maxIterations, double tolerance)
{
  double lambda0 = 0.0;
  double lambda1 = 0.0;
  double v1[n];
  double error = 10.0 * tolerance;
  int iter = 0;
  double norm;
  double Az[n * n];

  for (int i = 0; i < n; i++)
  {
    v1[i] = v0[i];
  }

  while (error > tolerance && iter < maxIterations)
  {
    matrixVectorProduct(n, n, A, v0, v1);
    norm = l2Norm(n, v1);
    for (int i = 0; i < n; i++)
    {
      v1[i] /= norm;
    }
    matrixVectorProduct(n, n, A, v1, Az);
    lambda1 = dotProduct(n, v1, Az);
    error = fabs(lambda1 - lambda0);
    lambda0 = lambda1;
    for (int i = 0; i < n; i++)
    {
      v0[i] = v1[i];
    }
    iter++;
  }
  return lambda0;
}
```

To test it, I randomly generated a diagonally dominant 12x12 matrix. I then used the power method to compute its largest eigenvalue and eigenvector. My test code is as follows:

```c
int main()
{
  int n = 12;
  double A[n * n];
  double v0[n];

  for (int i = 0; i < n; i++)
  {
    v0[i] = 1.0;
    for (int j = 0; j < n; j++)
    {
      if (i == j)
      {
        A[i * n + j] = 20 * n + rand() % 100 + 10;
      }
      else
      {

        A[i * n + j] = rand() % 20;
      }
    }
  }

  double lambdaMax = computePowerMethod(A, n, v0, 100, 0.0001);

  printf("A:\n");
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf("%8.2f ", A[i * n + j]);
    }
    printf("\n");
  }

  printf("\nv0:\n");
  for (int i = 0; i < n; i++)
  {
    printf("%f ", v0[i]);
  }
  printf("\n\nlambdaMax = %f\n", lambdaMax);
}
```

Output of test code:

```bash
A:
  257.00     9.00    13.00    18.00    10.00    12.00     4.00    18.00     3.00     9.00     0.00     5.00 
   12.00   292.00     7.00     3.00     7.00     9.00     0.00    12.00     3.00     9.00     9.00    17.00 
    0.00    13.00   349.00    18.00    16.00    15.00    17.00     6.00    12.00     7.00    10.00    13.00 
   19.00     9.00    19.00   271.00     7.00    12.00    13.00    16.00     5.00     5.00     8.00    11.00 
   14.00    17.00     1.00    13.00   258.00     4.00     8.00    10.00     4.00    16.00    10.00     3.00 
    2.00     6.00     9.00     4.00     1.00   303.00    17.00     8.00     8.00    13.00    18.00     1.00 
   15.00    13.00     5.00    14.00     3.00    16.00   275.00     9.00    15.00    14.00     9.00     1.00 
   17.00    15.00     5.00     4.00    11.00    18.00     8.00   273.00     5.00     2.00    12.00     6.00 
   16.00    17.00    18.00     4.00     1.00    17.00    11.00     8.00   287.00    18.00    17.00    17.00 
   14.00     4.00     9.00    11.00     5.00    15.00    15.00    18.00     2.00   349.00     8.00    12.00 
    0.00    17.00    14.00     8.00    15.00     8.00    13.00    10.00     6.00     2.00   292.00     5.00 
    2.00    12.00     7.00     1.00    15.00    12.00    11.00     1.00    10.00    11.00    18.00   307.00 

v0:
0.182998 0.200234 0.512750 0.241955 0.170041 0.238789 0.229235 0.179139 0.335177 0.451437 0.221033 0.263430

lambdaMax = 409.549041
```
