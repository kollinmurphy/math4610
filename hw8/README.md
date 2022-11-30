# Tasksheet 8 by Kollin Murphy

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
  double Az[n];

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

## Task 3: Reorganizing the algorithm for efficiency

I implemented the more efficient version of the power method that uses only one computation of the matrix-vector product per iteration. The code is as follows:

```c
double optimizedPowerMethod(double *A, int n, double *v0, int maxIterations, double tolerance)
{
  double lambda0 = 0.0;
  double lambda1 = 0.0;
  double y[n];
  double z[n];
  double w[n];
  double error = 10.0 * tolerance;
  int iter = 0;
  double norm;

  matrixVectorProduct(n, n, A, v0, y);

  while (error > tolerance && iter < maxIterations)
  {
    norm = l2Norm(n, y);
    for (int i = 0; i < n; i++)
    {
      z[i] = y[i] / norm;
    }
    matrixVectorProduct(n, n, A, z, w);
    lambda1 = dotProduct(n, z, w);
    error = fabs(lambda1 - lambda0);
    lambda0 = lambda1;
    for (int i = 0; i < n; i++)
    {
      y[i] = w[i];
    }
    iter++;
  }
  for (int i = 0; i < n; i++)
  {
    v0[i] = z[i];
  }
  return lambda0;
}
```

To compare the performance of the two versions of the power method, I generated random matrices of size 10x10, 100x100, and 200x200. I then computed the largest eigenvalue and eigenvector of each matrix using both versions of the power method. I repeated this process 100 times for each matrix size. I then computed the average time taken by each version of the power method to compute the largest eigenvalue and eigenvector of each matrix size. My test code is as follows:


```c
int main()
{
  int n = 10;
  double A[n * n];
  double v0[n];

  double timeRegular;
  double timeOptimized;
  double start;
  double end;

  int m = 100;

  double maxIter = 100;
  double tolerance = 0.0001;
  double result;

  printf("maxIter = %f; tolerance = %f\n\n", maxIter, tolerance);

  printf("10x10\n");

  timeRegular = 0;
  timeOptimized = 0;
  for (int i = 0; i < m; i++)
  {
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
    start = omp_get_wtime();
    result = computePowerMethod(A, n, v0, maxIter, tolerance);
    end = omp_get_wtime();
    timeRegular += end - start;

    start = omp_get_wtime();
    result = optimizedPowerMethod(A, n, v0, maxIter, tolerance);
    end = omp_get_wtime();
    timeOptimized += end - start;
  }
  printf("regular:   %f\n", timeRegular / m);
  printf("optimized: %f\n", timeOptimized / m);

  printf("\n100x100\n");
  n = 100;
  for (int i = 0; i < m; i++)
  {
    double A2[n * n];
    double v02[n];
    for (int i = 0; i < n; i++)
    {
      v02[i] = 1.0;
      for (int j = 0; j < n; j++)
      {
        if (i == j)
        {
          A2[i * n + j] = 20 * n + rand() % 100 + 10;
        }
        else
        {
          A2[i * n + j] = rand() % 20;
        }
      }
    }
    start = omp_get_wtime();
    result = computePowerMethod(A2, n, v02, maxIter, tolerance);
    end = omp_get_wtime();
    timeRegular += end - start;

    start = omp_get_wtime();
    result = optimizedPowerMethod(A2, n, v02, maxIter, tolerance);
    end = omp_get_wtime();
    timeOptimized += end - start;
  }
  printf("regular:   %f\n", timeRegular / m);
  printf("optimized: %f\n", timeOptimized / m);


  printf("\n200x200\n");
  n = 200;

  for (int i = 0; i < m; i++) {
    double A3[n * n];
    double v03[n];
    for (int i = 0; i < n; i++)
    {
      v03[i] = 1.0;
      for (int j = 0; j < n; j++)
      {
        if (i == j)
        {
          A3[i * n + j] = 20 * n + rand() % 100 + 10;
        }
        else
        {
          A3[i * n + j] = rand() % 20;
        }
      }
    }
    start = omp_get_wtime();
    result = computePowerMethod(A3, n, v03, maxIter, tolerance);
    end = omp_get_wtime();
    timeRegular += end - start;

    start = omp_get_wtime();
    result = optimizedPowerMethod(A3, n, v03, maxIter, tolerance);
    end = omp_get_wtime();
    timeOptimized += end - start;
  }
  printf("regular:   %f\n", timeRegular / m);
  printf("optimized: %f\n", timeOptimized / m);
}
```

The results of my test are as follows:

```
maxIter = 100.000000; tolerance = 0.000100

10x10
regular:   0.000080
optimized: 0.000005

100x100
regular:   0.003056
optimized: 0.001078

200x200
regular:   0.018446
optimized: 0.012424
```

As can be seen, the optimized version of the power method is significantly faster than the regular version. The optimized version of the power method is approximately 16 times faster than the regular version for a 10x10 matrix, approximately 3 times faster for a 100x100 matrix, and approximately 1.5 times faster for a 200x200 matrix.

The results also varied significantly based on the given tolerance. When the tolerance was larger, the optimized algorithm outperformed the other algorithm by even more on the larger matrices. Here is a result with a tolerance of 0.01:

```
maxIter = 100.000000; tolerance = 0.010000

10x10
regular:   0.000016
optimized: 0.000003

100x100
regular:   0.001447
optimized: 0.000222

200x200
regular:   0.007309
optimized: 0.001034
```

The optimized version of the power method is approximately 5 times faster than the regular version for a 10x10 matrix, approximately 6.5 times faster for a 100x100 matrix, and approximately 7 times faster for a 200x200 matrix.

## Task 4: Parallelizing the Power Method

I parallelized the power method by using OpenMP to parallelize the matrix-vector multiplication. I used the following code to parallelize the matrix-vector multiplication:

```c
void parallelizedMatrixVectorProduct(int n, int m, double *matrix, double *vector, double *result)
{
  double s;
#pragma omp parallel for
  for (int i = 0; i < n; i++)
  {
    result[i] = 0;
    for (int j = 0; j < m; j++)
    {
      result[i] += matrix[i * n + j] * vector[j];
    }
  }
}
```

I noticed that the parallelized version was faster on larger matrices, but slower on smaller matrices. This is to be expected due to the overhead of spinning up new threads when there aren't many operations to perform. Because of this, I also added code to test 500x500 matrices. Here are the results of my test:

```
maxIter = 100.000000; tolerance = 0.000100

10x10
optimized:    0.000032
parallelized: 0.000241

100x100
optimized:    0.001586
parallelized: 0.001679

200x200
optimized:    0.009567
parallelized: 0.007905

500x500
optimized:    0.087857
parallelized: 0.040290
```

Around the size of a 200x200 matrix, the parallelized version started to outperform the non-parallelized implementation. The parallelized version is approximately 7.5 times slower than the optimized version for a 10x10 matrix and approximately 1.3 times slower for a 100x100 matrix. The parallelized version is approximately 1.2 times faster for a 200x200 matrix and approximately 2.2 times faster for a 500x500 matrix.

## Task 5: Implementing Jacobi Iteration

I implemented Jacobi iteration to approximate the solution to a system of linear equations using the residual formula. I used the following code to implement Jacobi iteration:

```c
void jacobiIteration(int n, double *A, double *b, double *x0, double tolerance, int maxIterations)
{
  double error = 10.0 * tolerance;
  int iter = 0;
  double v[n];
  double r[n];
  double x1[n];

  for (int i = 0; i < n; i++)
    v[i] = x0[i];

  for (int i = 0; i < n; i++)
  {
    double sum = b[i];
    for (int j = 0; j < n; j++)
      sum -= A[i * n + j] * v[j];
    r[i] = sum;
  }

  while (error > tolerance && iter < maxIterations)
  {
    // get next x
    for (int i = 0; i < n; i++)
      x1[i] = v[i] + r[i] / A[i * n + i];

    // get next error
    error = 0;
    for (int i = 0; i < n; i++)
    {
      double diff = x1[i] - v[i];
      error += diff * diff;
    }
    error = sqrtf(error);

    iter++;

    // get next residual
    for (int i = 0; i < n; i++)
    {
      double sum = b[i];
      for (int j = 0; j < n; j++)
        sum -= A[i * n + j] * x1[j];
      r[i] = sum;
    }

    // overwrite v
    for (int i = 0; i < n; i++)
      v[i] = x1[i];
  }

  for (int i = 0; i < n; i++)
    x0[i] = v[i];
}
```

I generated a random 100x100, diagonally dominant matrix and a random 100x1 vector. I then used the Jacobi iteration method to approximate the solution to the system of linear equations. I used the following code to test the Jacobi iteration method:

```c

int main()
{
  int n = 100;
  double A[n * n];
  double b[n];
  double x0[n];
  for (int i = 0; i < n; i++)
  {
    x0[i] = 0.0;
    b[i] = rand() % 10000;
    for (int j = 0; j < n; j++)
    {
      if (i == j)
        A[i * n + j] = 21 * n + rand() % 100 + 10;
      else
        A[i * n + j] = rand() % 20 + 1;
    }
  }

  jacobiIteration(n, A, b, x0, 0.0001, 100);

  for (int i = 0; i < n; i++)
    printf("%f ", x0[i]);
  printf("\n");

  return 0;
}
```

The results of my test were as follows:

```
-0.007831 -0.006950 -0.006786 -0.006073 -0.005442 -0.005682 -0.005095 -0.004575 -0.004012 -0.003835 -0.001699 -0.003030 -0.003004 -0.000661 -0.001475 -0.001031 -0.000293 0.000181 -0.000122 0.000458 0.001698 0.001496 0.002705 0.002652 0.004035 0.005046 0.005114 0.004733 0.006122 0.006940 0.008057 0.006366 0.006763 0.008477 0.007860 0.007982 0.008367 0.010303 0.009504 0.011262 0.011392 0.011695 0.012126 0.011251 0.013569 0.012722 0.013204 0.014717 0.014678 0.014799 0.015722 0.016550 0.016127 0.016411 0.018608 0.017616 0.018390 0.019346 0.020792 0.020193 0.020634 0.021541 0.021843 0.021528 0.021599 0.022597 0.023053 0.023219 0.025099 0.024160 0.025379 0.025391 0.025677 0.026880 0.026755 0.028543 0.027957 0.028269 0.028360 0.028743 0.030025 0.029855 0.030727 0.032330 0.032343 0.031589 0.032944 0.033141 0.033929 0.032301 0.033539 0.034407 0.034214 0.035230 0.036885 0.038232 0.036557 0.036152 0.037461 0.038417
```

It was able to successfully approximate the solution to the system of linear equations. I also tested it on a smaller matrix with a known solution and it was able to successfully approximate the solution to the system of linear equations within the given tolerance.
