from math import e
from rootFinding import fixedPoint

def f(x):
  return x * pow(e, -x)

def g(x):
  return x - f(x)

print(fixedPoint(g, 0.5, 0.0000001, 1000))
