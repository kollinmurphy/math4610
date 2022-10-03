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
