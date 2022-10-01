#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fixedPoint(double (*g)(double), double initialApproximation, double tolerance, int maxIterations) {
  double error = 10.0 * tolerance;
  double x0 = initialApproximation;
  double x1 = 0;
  int iterations = 0;

  while (error > tolerance && iterations < maxIterations) {
    x1 = g(x0);
    error = fabs(x1 - x0);
    x0 = x1;
    iterations += 1;
  }

  return x0;
}
