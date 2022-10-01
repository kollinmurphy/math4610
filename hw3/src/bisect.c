#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double bisect(double (*f)(double), double a, double b, double tolerance) {
  double fa = f(a);
  double fb = f(b);

  if (fa * fb >= 0) {
    printf("ERROR: Intermediate Value Theorem is not satisfied by the bisection endpoints.");
    exit(1);
  }

  double k = ceil(-log(tolerance / (b - a)) / log(2));
  double c = a;
  double fc;

  for (int i = 0; i < k; i++) {
    c = (a + b) / 2;
    fc = f(c);
    if (fa * fc < 0) {
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }
  }
  return c;
}
