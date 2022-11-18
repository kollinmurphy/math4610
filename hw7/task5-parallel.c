#include <stdio.h>
#include <omp.h>

void outerProductOfVectors(int n, int m, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    for (int i = id; i < n; i = i + numThreads) {
      for (int j = 0; j < m; ++j) {
        c[i * m + j] = a[i] * b[j];
      }
    }
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}
