from math import pi, pow, e, log
import matplotlib.pyplot as plt

def simpsonsRule(f, a, b, n):
  h = (b - a) / n
  x = a
  sum = 0
  for i in range(n):
    sum += f(x) + 4 * f(x + h / 2) + f(x + h)
    x += h
  return sum * h / 6


def calculateNFromH(a, b, h):
  return (b - a) / h

hvals = [1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4196]
nvals = [int(calculateNFromH(0, 1, h)) for h in hvals]

# I calculated the actual value of the integral using Wolfram Alpha
actual = 0.6498803300786573037276521829129935239240253152760926685227030955

def f(x):
  return pow(e, -x * x)

approximations = [simpsonsRule(f, 0, pi / 4, n) for n in nvals]
errors = [abs(actual - approximation) for approximation in approximations]

plt.plot([log(h) for h in hvals], [log(error) for error in errors])
plt.title("Log-Log Plot of Error vs. h")
plt.xlabel("log(h)")
plt.ylabel("log(error)")
plt.savefig("task5.png")
plt.show()
