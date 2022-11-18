#include <stdio.h>
#include <omp.h>

void hadamardProductOfMatrices(int m, int n, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    for (int i = id; i < m; i = i + numThreads) {
      for (int j = 0; j < n; ++j) {
        c[i * n + j] = a[i * n + j] * b[i * n + j];
      }
    }
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}

int main() {
  int m = 100;
  int n = 1000;
  double a[m][n];
  double b[m][n];
  double c[m][n];
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i][j] = i + j + 1;
      b[i][j] = i + j + 51;
    }
  }
  printf("m = %d, n = %d\n", m, n);
  hadamardProductOfMatrices(m, n, &a[0][0], &b[0][0], &c[0][0]);
  return 0;
}
