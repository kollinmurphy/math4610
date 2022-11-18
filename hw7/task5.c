#include <stdio.h>
#include <omp.h>

void outerProductOfVectors(int n, int m, double *a, double *b, double *c) {
  double start = omp_get_wtime();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; ++j) {
      c[i * m + j] = a[i] * b[j];
    }
  }
  double end = omp_get_wtime();
  printf("Time: %.16f\n", end - start);
}

int main() {
  int n = 5;
  int m = 10;
  double a[n];
  double b[m];
  double c[n][m];
  for (int i = 0; i < n; ++i) {
    a[i] = i + 1;
  }
  for (int j = 0; j < m; ++j) {
    b[j] = j + 51;
  }
  printf("n = %d, m = %d\n", n, m);
  outerProductOfVectors(n, m, &a[0], &b[0], &c[0][0]);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      printf("%f ", c[i][j]);
    }
    printf("\n");
  }
  return 0;
}
