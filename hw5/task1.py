# setting up importing code from other homework directories
import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from hw2.newtons import newtons
import matplotlib.pyplot as plt

def implicitEulerMethod(alpha, beta, p0, n=100, h=1.0):
  """Solve the logistic equation using the implicit Euler method."""

  def fPrime(p):
    return 1 - alpha * h + 2 * h * beta * p

  ptable = [0 for i in range(n)]
  ptable[0] = p0
  for i in range(1,n):
    def f(p):
      return p - ptable[i - 1] - h * (alpha * p - beta * p * p)
    pval = newtons(f, fPrime, ptable[i - 1], 100, 1e-6)
    ptable[i] = pval
  return zip(*[[i*h, ptable[i]] for i in range(n)])
  
tests = [
  [0.2, 0.0005, 10.0, 75, 1.0, "Test 1", "task1-test1.png"],
  [0.01, 0.0005, 10.0, 350, 2.0, "Test 2", "task1-test2.png"],
  [2.0, 0.0005, 10.0, 75, 0.1, "Test 3", "task1-test3.png"],
]

for test in tests:
  tvals, pvals = implicitEulerMethod(test[0], test[1], test[2], test[3], test[4])
  plt.title(test[5])
  plt.xlabel('t')
  plt.ylabel('P(t)')
  plt.plot(tvals, pvals)
  plt.savefig(test[6])
  plt.clf()
