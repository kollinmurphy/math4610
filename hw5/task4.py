from math import e, pow, pi

def trapezoidRule(f, a, b, n):
  h = (b - a) / n
  sum = 0.5 * (f(a) + f(b))
  for i in range(1, n):
    sum += f(a + i * h)
  return sum * h

tests = [2, 4, 8, 16, 1000]

def f(x):
  return pow(e, -x * x)

for n in tests:
  print(f"n={n}; approx={trapezoidRule(f, 0, pi / 4, n)}")
