#include <stdio.h>
#include <math.h>

double f(double x) {
  return 1 / sqrtf(1 - (x * x));
}

double trapezoidRule(double (*f)(double), double a, double b, int n) {
  double h = (b - a) / n;
  double sum = 0.5 * (f(a) + f(b));
  for (int i = 1; i < n; ++i) {
    sum += f(a + i * h);
  }
  return sum * h;
}

int main() {
  double pi = 6 * trapezoidRule(f, 0, 0.5, 100000000);
  printf("pi = %.16f\n", pi);
  return 0;
}
