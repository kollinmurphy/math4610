#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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
