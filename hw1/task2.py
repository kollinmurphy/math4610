from math import e
from rootFinder import fixedPointRootFinder

def f(x):
  return x * pow(e, -x)

print(fixedPointRootFinder(f, 0.5, 0.0000001, 1000))
