#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double l2Norm(int n, double vector[n])
{
  double sum = 0;
  for (int i = 0; i < n; i++)
  {
    sum += vector[i] * vector[i];
  }
  return sqrtf(sum);
}

void jacobiIteration(int n, double *A, double *b, double *x0, double *x, double tolerance, int maxIterations)
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
    x[i] = v[i];
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

  jacobiIteration(n, A, x0, initialForJacobi, y, tolerance, maxIterations); // y = A^-1 * x0

  while (error > tolerance && iter < maxIterations)
  {
    norm = l2Norm(n, y);
    for (int i = 0; i < n; i++)
      z[i] = y[i] / norm;

    jacobiIteration(n, A, z, initialForJacobi, w, tolerance, maxIterations); // w = A^-1 * z
    lambda1 = dotProduct(n, z, w);                                           // lambda1 = z^T * w
    error = fabs(lambda1 - lambda0);
    lambda0 = lambda1;
    iter++;
    for (int i = 0; i < n; i++)
      y[i] = w[i];
  }

  for (int i = 0; i < n; i++)
    result[i] = y[i];
  return 1.0 / lambda1;
}

int main()
{
  int n = 3;
  double A[9] = {67, 8, 4, 8, 90, 8, 4, 8, 83};
  double v[n];
  double x0[n];
  for (int i = 0; i < n; i++)
  {
    x0[i] = 1.0;
  }

  double mu = 77.0;
  for (int i = 0; i < n; i++)
    A[i * n + i] -= mu;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%f ", A[i * n + j]);
    }
    printf("\n");
  }
  printf("\n\n");
  double lambda = inversePowerMethod(n, A, v, x0, 1e-8, 100) + mu;
  printf("lambda = %f\n", lambda);
  printf("v = [");
  for (int i = 0; i < n; i++)
  {
    printf("%f", v[i]);
    if (i < n - 1)
      printf(", ");
  }
  printf("]\n");
  return 0;
}