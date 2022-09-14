from sys import argv
from math import ceil, log

verbose = False
for arg in argv:
  if arg == '-v' or arg == '--v' or arg == '--verbose':
    verbose = True

def fixedPoint(g, initial_approximation, tolerance, max_iterations):
  error = 10.0 * tolerance
  x0 = initial_approximation
  x1 = 0
  iterations = 0
  
  if verbose:
    print("{:<12} {:<24} {:<24}".format("Iteration", "Approx. Root", "Abs. Error"))

  while error > tolerance and iterations < max_iterations:
    x1 = g(x0)
    error = abs(x1 - x0)
    x0 = x1
    iterations += 1

    if verbose:
      print("{:<12} {:<24} {:<24}".format(iterations, x0, error))
  return x0


def bisect(f, a, b, tolerance):
  fa = f(a)
  fb = f(b)

  if fa * fb >= 0:
    raise Exception('Intermediate Value Theorem is not satisfied by the initial conditions')

  k = ceil(-log(tolerance / (b - a)) / log(2))
  c = a

  if verbose:
    print("{:<12} {:<24} {:<24}".format("Iteration", "Approx. Root", "Abs. Error"))

  for i in range(0, k):
    c = (a + b) / 2
    fc = f(c)
    if fa * fc < 0:
      b = c
      fb = fc
    else:
      a = c
      fa = fc
    if verbose:
      error = abs(a - b)
      print("{:<12} {:<24} {:<24}".format(i, c, error))
  
  return c
