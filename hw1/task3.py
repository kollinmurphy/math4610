from math import cos, e, pi
from rootFinding import fixedPoint

def f(x):
  return 10.14 * pow(e, x * x) * cos(pi / x)

def g(x):
  return x - f(x)

for i in range(-12, 29):
  print()
  print(i/4)
  try:
    print(fixedPoint(g, i/4, 0.01, 100))
  except:
    print('failure')
