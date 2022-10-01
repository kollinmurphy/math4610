#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int HYBRID_NEWTONS_BISECTION_ITERATIONS = 4;

double hybridNewtons(double (*f)(double), double (*fPrime)(double), double a, double b, double tolerance, int maxIterations) {
  double error = 10.0 * tolerance;
  int i = 0;

  double x0 = (a + b) / 2;
  double x1, e1, fa, fb, fc, c;
  while (error > tolerance && i < maxIterations) {
    x1 = x0 - (f(x0) / fPrime(x0));
    e1 = fabs(x1 - x0);
    
    if (e1 > error) {
      fa = f(a);
      for (int i = 0; i < HYBRID_NEWTONS_BISECTION_ITERATIONS; i++) {
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
      x0 = (a + b) / 2;
      error = fabs(a -b);
    } else {
      x0 = x1;
      error = e1;
    }
    i++;
  }
  return x0;
}