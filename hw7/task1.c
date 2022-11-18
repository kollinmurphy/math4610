#include <stdio.h>
#include <omp.h>

void matrixProduct(int aRows, int aCols, int bCols, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  #pragma omp parallel for
  for (int i = 0; i < aRows; ++i) {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    printf("Thread %d of %d working on iteration %d\n", id, numThreads, i);
    for (int j = 0; j < bCols; ++j) {
      double sum = 0;
      for (int k = 0; k < aCols; ++k) {
        sum += a[i * aCols + k] * b[k * bCols + j];
      }
      c[i * bCols + j] = sum;
    }
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}

int main() {
  double a[2][3] = {
    {1, 2, 3},
    {4, 5, 6}
  };
  double b[3][2] = {
    {7, 8},
    {9, 10},
    {11, 12}
  };
  double c[2][2];
  matrixProduct(2, 3, 2, &a[0][0], &b[0][0], &c[0][0]);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      printf("%f ", c[i][j]);
    }
    printf("\n");
  }
  return 0;
}
