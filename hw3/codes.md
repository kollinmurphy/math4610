# C Codes

[Return to tasksheet answers](README.md)

## Bisect

Implementation of `bisect.c`.

```c
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
```

## Fixed Point

Implementation of `fixed_point.c`.

```c
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
```

## Hybrid-Newtons

Implementation of `hybrid_newtons.c`.

```c
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
```

## Hybrid-Secant

Implementation of `hybrid_secant.c`.

```c
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
```

## Newtons

Implementation of `newtons.c`.

```c
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
```

## Secant

Implementation of `secant.c`.

```c
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
```

## Test Code

I tested my various functions with the following code in `test_root_finding.c`.

```c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x) {
  return 10.0 * sin(x);
}

double xex(double x) {
  return x * pow(M_E, -x);
}

double xexPrime(double x) {
  return pow(M_E, -x) - (x * pow(M_E, -x));
}

double g(double x) {
  return x - (x * pow(M_E, -x));
}

double fPrime(double x) {
  return 10.0 * cos(x);
}

double tolerance = 0.0001;
int maxIterations = 200;

double fixedPoint(double (*g)(double), double initialApproximation, double tolerance, int maxIterations);
double bisect(double (*f)(double), double a, double b, double tolerance);
double hybridNewtons(double (*f)(double), double (*fPrime)(double), double a, double b, double tolerance, int maxIterations);
double hybridSecant(double (*f)(double), double x0, double x1, double tolerance, int maxIterations);
double newtons(double (*f)(double), double (*fPrime)(double), double x, double tolerance, int maxIterations);
double secant(double (*f)(double), double x0, double x1, double tolerance, int maxIterations);

int main() {
  printf("fixedPoint\n\tactual:\t\t%.32f\n\texpected:\t%.32f\n\n", fixedPoint(g, 0.5, tolerance, maxIterations), 0.0);
  printf("bisect\n\tactual:\t\t%.32f\n\texpected:\t%.32f\n\n", bisect(f, 2.0, 4.0, tolerance), M_PI);
  printf("newtons\n\tactual:\t\t%.32f\n\texpected:\t%.32f\n\n", newtons(xex, xexPrime, -1.0, tolerance, maxIterations), M_PI);
  printf("secant\n\tactual:\t\t%.32f\n\texpected:\t%.32f\n\n", secant(f, 2.0, 3.0, tolerance, maxIterations), M_PI);
  printf("hybridNewtons\n\tactual:\t\t%.32f\n\texpected:\t%.32f\n\n", hybridNewtons(f, fPrime, 2.0, 4.0, tolerance, maxIterations), M_PI);
  printf("hybridSecant\n\tactual:\t\t%.32f\n\texpected:\t%.32f\n\n", hybridSecant(f, 2.0, 4.0, tolerance, maxIterations), M_PI);
}
```
