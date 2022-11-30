#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
      {
        A[i * n + j] = 21 * n + rand() % 100 + 10;
      }
      else
      {
        A[i * n + j] = rand() % 20 + 1;
      }
    }
  }

  jacobiIteration(n, A, b, x0, 0.0001, 100);

  for (int i = 0; i < n; i++)
  {
    printf("%f ", x0[i]);
  }
  printf("\n");

  return 0;
}
