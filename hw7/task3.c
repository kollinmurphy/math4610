#include <stdio.h>
#include <omp.h>

void hadamardProduct(int n, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    for (int i = id; i < n; i = i + numThreads) {
      c[i] = a[i] * b[i];
    }
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}

int main() {
  double a[100];
  double b[100];
  double c[100];
  for (int i = 0; i < 100; ++i) {
    a[i] = i + 1;
    b[i] = i + 51;
  }
  hadamardProduct(100, &a[0], &b[0], &c[0]);
  for (int i = 0; i < 100; ++i) {
    printf("%f ", c[i]);
  }
  printf("\n");
  return 0;
}
