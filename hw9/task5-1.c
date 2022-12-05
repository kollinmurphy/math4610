#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

double l2Norm(int n, double vector[n])
{
  double sum = 0;
  for (int i = 0; i < n; i++)
  {
    sum += vector[i] * vector[i];
  }
  return sqrtf(sum);
}

void gaussSeidel(int n, double *A, double *b, double *initialApproximation, double *result, double tolerance, int maxIterations)
{
  double x0[n];
  for (int i = 0; i < n; i++)
    x0[i] = initialApproximation[i];

  for (int iter = 0; iter < maxIterations; iter++)
  {
    double x1[n];
    for (int i = 0; i < n; i++)
      x1[i] = 0;

    for (int i = 0; i < n; i++)
    {
      double s1 = 0.0;
      for (int j = 0; j < i; j++)
        s1 += A[i * n + j] * x1[j];
      double s2 = 0.0;
      for (int j = i + 1; j < 4; j++)
        s2 += A[i * n + j] * x0[j];
      x1[i] = (b[i] - s1 - s2) / A[i * n + i];
    }
    int converged = 1;
    for (int i = 0; i < n; i++)
    {
      if (fabs(x1[i] - x0[i]) > tolerance)
      {
        converged = 0;
        break;
      }
    }
    for (int i = 0; i < n; i++)
      x0[i] = x1[i];
    if (converged == 1)
      break;
  }

  double error[n];
  for (int i = 0; i < n; i++)
  {
    error[i] = 0.0;
    for (int j = 0; j < n; j++)
    {
      error[i] += A[i * n + j] * x0[j];
    }
    error[i] -= b[i];
  }

  for (int i = 0; i < n; i++)
    result[i] = x0[i];
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

void matrixVectorProduct(int n, int m, double *matrix, double *vector, double *result)
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

double inversePowerMethod(int n, double *A, double *result, double *x0, double tolerance, int maxIterations)
{
  double lambda0 = 0.0;
  double lambda1 = 0.0;
  double error = 10.0 * tolerance;
  int iter = 0;
  double y[n];
  double z[n];
  double w[n];
  double initialForJacobi[n];
  double norm;

  for (int i = 0; i < n; i++)
  {
    initialForJacobi[i] = 1.0;
    y[i] = 0.0;
    w[i] = 0.0;
  }

  gaussSeidel(n, A, x0, initialForJacobi, y, tolerance, maxIterations); // y = A^-1 * x0

  while (error > tolerance && iter < maxIterations)
  {
    norm = l2Norm(n, y);
    for (int i = 0; i < n; i++)
      z[i] = y[i] / norm;

    gaussSeidel(n, A, z, initialForJacobi, w, tolerance, maxIterations); // w = A^-1 * z
    lambda1 = dotProduct(n, z, w);                                           // lambda1 = z^T * w
    error = fabs(lambda1 - lambda0);
    lambda0 = lambda1;
    iter++;
    for (int i = 0; i < n; i++)
      y[i] = w[i];
  }

  if (iter == maxIterations)
  {
    printf("Inverse power method did not converge in %d iterations\n", iter);
  }

  for (int i = 0; i < n; i++)
    result[i] = y[i];
  return 1.0 / lambda1;
}

double powerMethod(int n, double *A, double *result, double *v0, double tolerance, int maxIterations)
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
    result[i] = z[i];
  }
  return lambda0;
}

void copyMatrix(int n, double *A, double *B)
{
  for (int i = 0; i < n * n; i++)
    B[i] = A[i];
}

int main()
{
  int n = 3;
  double A[n * n];
  double B[n * n];
  double v[n];
  double x0[n];

  for (int i = 0; i < n; i++)
  {
    x0[i] = 1.0;
    for (int j = 0; j < n; j++)
    {
      if (i == j)
        A[i * n + j] = 20 * n + (rand() % 100);
      else {
        A[i * n + j] = rand() % 10;
        A[j * n + i] = A[i * n + j];
      }
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%f ", A[i * n + j]);
    }
    printf("\n");
  }
  printf("\n\n");

  double tolerance = 0.0001;
  int maxIterations = 100;

  double lambdaMax = powerMethod(n, A, v, x0, tolerance, maxIterations);
  printf("lambdaMax = %f\n", lambdaMax);

  double lambdaMin = inversePowerMethod(n, A, v, x0, tolerance, maxIterations);
  printf("lambdaMin = %f\n\n", lambdaMin);

  double lambdaRange = lambdaMax - lambdaMin;
  int partitions = n;
  double partitionSize = lambdaRange / partitions;

  #pragma omp parallel
  {
    int threadId = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    for (int i = 1 + threadId; i < partitions - 1; i += numThreads)
    {
      copyMatrix(n, A, B);
      double mu = lambdaMin + i * partitionSize;
      for (int i = 0; i < n; i++)
        B[i * n + i] -= mu;
      double lambda = inversePowerMethod(n, B, v, x0, tolerance, maxIterations);
      printf("lambda = %f\n", lambda + mu);
    }
  }

  return 0;
}
