from math import exp

def analyticSolution(alpha, beta, p0, n=100, h=1.0):
  """Solve the logistic equation using the analytic solution."""
  c = alpha / p0 - beta
  return zip(*[[i*h, alpha / (c * exp(-alpha * i*h) + beta)] for i in range(n)])
