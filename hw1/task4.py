from rootFinding import bisect

def f(x):
  return x * x - 4.1

print(bisect(f, a=-1, b=3, tolerance=0.000001))
