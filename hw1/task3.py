from math import cos, e, pi
from rootFinder import fixedPointRootFinder

def f(x):
  return 10.14 * pow(e, x * x) * cos(pi / x)

for i in range(-12, 29):
  print()
  print(i/4)
  try:
    print(fixedPointRootFinder(f, i/4, 0.01, 100))
  except:
    print('failure')
