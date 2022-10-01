#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double secant(double (*f)(double), double x0, double x1, double tolerance, int maxIterations) {
  int i = 0;
  double f0 = f(x0);
  double f1 = f(x1);
  double x2;

  while (fabs(f1) > tolerance && i < maxIterations) {
    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    x0 = x1;
    x1 = x2;
    f0 = f1;
    f1 = f(x1);
    i++;
  }

  return x1;
}
