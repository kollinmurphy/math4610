from math import e, pow, pi

def trapezoidRule(f, a, b, n):
  h = (b - a) / n
  sum = 0.5 * (f(a) + f(b))
  for i in range(1, n):
    sum += f(a + i * h)
  return sum * h

tests = [2,4,8,16, 100000]

def f(x):
  return pow(e, -x * x)

for test in tests:
  print(f"n={test}; approx={trapezoidRule(f, 0, pi / 4, test)}")
