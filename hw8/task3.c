#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

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
  double tolerance = 0.0000001;
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
