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
