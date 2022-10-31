# setting up importing code from other homework directories
import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from hw4.task5 import explicitEulerLogistic
from task1 import implicitEulerMethod
from task2 import analyticSolution
from matplotlib import pyplot as plt

tests = [
  [0.2, 0.0005, 10.0, 50, 2, "Test 1 (n=50, h=2.0)", "task3-test1-a.png"],
  [0.2, 0.0005, 10.0, 200, 0.5, "Test 1 (n=200, h=0.5)", "task3-test1-b.png"],
  [0.01, 0.0005, 10.0, 50, 15.0, "Test 2 (n=50, h=15.0)", "task3-test2-a.png"],
  [0.01, 0.0005, 10.0, 200, 3.75, "Test 2 (n=200, h=3.75)", "task3-test2-b.png"],
  [2.0, 0.0005, 10.0, 50, 0.2, "Test 3 (n=50, n=0.2)", "task3-test3-a.png"],
  [2.0, 0.0005, 10.0, 200, 0.05, "Test 3 (n=200, 0.05)", "task3-test3-b.png"],
]

for test in tests:
  tImplicit, pImplicit = implicitEulerMethod(test[0], test[1], test[2], test[3], test[4])
  tExplicit, pExplicit = explicitEulerLogistic(test[0], test[1], test[2], test[3], test[4])
  tAnalytic, pAanalytic = analyticSolution(test[0], test[1], test[2], test[3], test[4])

  plt.title(test[5])
  plt.xlabel('t')
  plt.ylabel('P(t)')

  plt.plot(tImplicit, pImplicit, 'r', label='Implicit Euler')
  plt.plot(tExplicit, pExplicit, 'g', label='Explicit Euler')
  plt.plot(tAnalytic, pAanalytic, 'b', label='Analytic')

  plt.legend(loc="lower right")
  plt.savefig(test[6])
  plt.clf()
