from math import cos, e, pi, sin
from sys import argv

verbose = False
for arg in argv:
  if arg == '-v' or arg == '--v' or arg == '--verbose':
    verbose = True

bisection_iterations = 4

def hybrid_newtons(f, fPrime, a, b, maxIterations, tolerance):
  error = 10.0 * tolerance
  i = 0

  if verbose:
    print("{:<12} {:<12} {:<24} {:<24}".format("Iteration", "Algorithm", "Approx. Root", "Abs. Error"))

  x0 = (a + b) / 2
  while (error > tolerance and i < maxIterations):
    x1 = x0 - (f(x0) / fPrime(x0))
    e1 = abs(x1 - x0)

    if (e1 > error):
      fa = f(a)
      for _ in range(bisection_iterations):
        c = (a + b) / 2
        fc = f(c)
        if fa * fc < 0:
          b = c
          fb = fc
        else:
          a = c
          fa = fc
      x0 = (a + b) / 2
      error = abs(a - b)
      if verbose:
        print("{:<12} {:<12} {:<24} {:<24}".format(i, 'bisection', x0, error))
    else:
      x0 = x1
      error = e1
      if verbose:
        print("{:<12} {:<12} {:<24} {:<24}".format(i, 'newtons', x0, error))
    i += 1
  return x0

def f(x):
  return 10.14 * pow(e, pow(x, 2)) * cos(pi / x)

def fPrime(x):
  return ((1014 * cos(pi / x) * pow(x, 3)) + (507 * pi * sin(pi / x))) * pow(e, pow(x, 2)) / (50 * pow(x, 2))

print(hybrid_newtons(f, fPrime, -3, 7, 10, 0.001))

# def f(x):
#   return 1000 * sin(x)

# def fPrime(x):
#   return 1000 * cos(x)

# print(hybrid_newtons(f, fPrime, -2.1, 1, 100, 10e-12))
