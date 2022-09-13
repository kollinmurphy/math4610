
from math import e


def fixedPointRootFinder(f, initial_approximation, tolerance, max_iterations):
  error = 10.0 * tolerance
  x0 = initial_approximation
  x1 = 0
  iterations = 0

  def g(x):
    return x - f(x)

  while error > tolerance and iterations < max_iterations:
    x1 = g(x0)
    error = abs(x1 - x0)
    x0 = x1
    iterations += 1
  return x0

def f(x):
  return x * pow(e, -x)

print(fixedPointRootFinder(f, 0.5, 0.0000001, 1000))
