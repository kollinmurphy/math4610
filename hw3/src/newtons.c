#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double newtons(double (*f)(double), double (*fPrime)(double), double x, double tolerance, int maxIterations) {
  int i = 0;
  double y = 10.0 * tolerance;

  while (fabs(y) > tolerance && i < maxIterations) {
    y = f(x) / fPrime(x);
    x = x - y;
    i++;
  }

  return x;
}
