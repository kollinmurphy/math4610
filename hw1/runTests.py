from rootFinding import fixedPoint, bisect

def f(x):
  return x * x - 4.1

def g(x):
  return x - 0.1 * f(x)

tolerance = 0.000001
max_iterations = 10000

print("Finding roots of f(x) = x^2 - 4.1")
print()
print(f"Fixed point: {fixedPoint(g, initial_approximation=-1, tolerance=tolerance, max_iterations=max_iterations)}")
print(f"Bisection: {bisect(f, a=-1, b=3, tolerance=tolerance)}")
