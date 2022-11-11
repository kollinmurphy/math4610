#include <stdio.h>
#include <math.h>
#include <omp.h>

double f(double x) {
  return 1 / sqrtf(1 - (x * x));
}

static long num_steps = 100000000;
#define NUM_THREADS 4

double trapezoidRule(double (*f)(double), double a, double b, int n) {
  double h = (b - a) / n;
  double sum[NUM_THREADS];
  int actualNumThreads = 0;
  double accumulator = 0.5 * (f(a) + f(b));

  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    if (id == 0) actualNumThreads = numThreads;

    for (int i = 1 + id; i < n; i = i + numThreads) {
      sum[id] += f(a + i * h);
    }
  }

  for (int i = 0; i < actualNumThreads; ++i) {
    accumulator += sum[i];
  }
  return accumulator * h;
}

int main() {
  double pi = 6 * trapezoidRule(f, 0, 0.5, num_steps);
  printf("pi = %.16f\n", pi);
  return 0;
}
