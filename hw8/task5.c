#include <math.h>

void jacobiIteration(int n, double *A, double *b, double *x0, double tolerance, int maxIterations)
{
  double error = 10.0 * tolerance;
  int iter = 0;
  double v[n];
  double r[n];
  double x1[n];

  for (int i = 0; i < n; i++)
  {
    v[i] = x0[i];
  }

  for (int i = 0; i < n; i++)
  {
    double sum = b[i];
    for (int j = 0; j < n; j++)
    {
      // if (i != j) {
      sum -= A[i * n + j] * v[j]; // or maybe v[i]?
      // }
    }
    r[i] = sum;
  }

  while (error > tolerance && iter < maxIterations)
  {
    // get next x
    for (int i = 0; i < n; i++)
    {
      x1[i] = v[i] + r[i] / A[i * n + i];
    }

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
      {
        // if (i != j) {
        sum -= A[i * n + j] * x1[j]; // or maybe x1[i]?
        // }
      }
      r[i] = sum;
    }

    // overwrite v
    for (int i = 0; i < n; i++)
    {
      v[i] = x1[i];
    }
  }
}
