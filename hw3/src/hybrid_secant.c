#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int HYBRID_SECANT_BISECTION_ITERATIONS = 4;

double hybridSecant(double (*f)(double), double x0, double x1, double tolerance, int maxIterations) {
  double lastA = x0;
  double lastB = x1;
  double error = 10.0 * tolerance;
  int i = 0;
  double f0, f1, x2, e, a, b, c, fa, fb, fc;

  while (error > tolerance && i < maxIterations) {
    f0 = f(x0);
    f1 = f(x1);
    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    e = fabs(x2 - x1);

    if (e > error) {
      a = x0;
      b = x1;
      fa = f(a);
      fb = f(b);
      if (fa * fb > 0) {
        a = lastA;
        b = lastB;
        fa = f(a);
        fb = f(b);
        if (fa * fb > 0) {
          printf("ERROR: Intermediate Value Theorem is not satisfied by the bisection endpoints.");
          exit(1);
        }
      }
      for (int i = 0; i < HYBRID_SECANT_BISECTION_ITERATIONS; i++) {
        c = (a + b) / 2;
        fc = f(c);
        if (fa * fc < 0) {
          b = c;
          fb = fc;
        } else {
          a = c;
          fa = fc;
        }
        x0 = x1;
        x1 = (a + b) / 2;
        lastA = a;
        lastB = b;
      }
    } else {
      x0 = x1;
      x1 = x2;
      error = e;
    }
    i++;
  }
  return x1;
}
