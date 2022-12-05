#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

    printf("Iteration %d: [ ", iter);
    for (int i = 0; i < n; i++)
      printf("%f ", x0[i]);
    printf("]\n");

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

  printf("Solution: [ ");
  for (int i = 0; i < n; i++)
    printf("%f ", x0[i]);
  printf("]\n");

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

  printf("Error: [ ");
  for (int i = 0; i < n; i++)
    printf("%f ", error[i]);
  printf("]\n");

  for (int i = 0; i < n; i++)
    result[i] = x0[i];
}

int main()
{
  double tolerance = 0.00001;
  int maxIterations = 1000;
  int n = 4;
  double A[16] = {
      10.0, -1.0, 2.0, 0.0,
      -1.0, 11.0, -1.0, 3.0,
      2.0, -1.0, 10.0, -1.0,
      0.0, 3.0, -1.0, 8.0};
  double b[4] = {6.0, 25.0, -11.0, 15.0};
  double x0[4] = {0.0, 0.0, 0.0, 0.0};
  double result[4];
  gaussSeidel(n, A, b, x0, result, tolerance, maxIterations);
  return 0;
}
