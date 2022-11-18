#include <stdio.h>
#include <omp.h>

void hadamardProduct(int n, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  for (int i = 0; i < n; ++i) {
    c[i] = a[i] * b[i];
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}

int main() {
  double a[4] = {1, 2, 3, 4};
  double b[4] = {5, 6, 7, 8};
  double c[4];
  hadamardProduct(4, &a[0], &b[0], &c[0]);
  for (int i = 0; i < 4; ++i) {
    printf("%f ", c[i]);
  }
  printf("\n");
  return 0;
}
