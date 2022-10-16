from math import pi, tan

def approxSecondDerivative(f, x0, h):
  return f(x0 - h) - 2 * f(x0) + f(x0 + h)

def testFunction(x):
  return ((x - pi / 2) * tan(x) * tan(x)) / (x * x + 65)

print(f"h = 1 \t\tapprox = \t {approxSecondDerivative(testFunction, pi / 4, 1)}")
print(f"h = 0.1 \tapprox =\t{approxSecondDerivative(testFunction, pi / 4, 0.1)}")
print(f"h = 0.001 \tapprox = \t{approxSecondDerivative(testFunction, pi / 4, 0.001)}")
print(f"h = 0.00001 \tapprox = \t{approxSecondDerivative(testFunction, pi / 4, 0.00001)}")
print(f"h = 0.0000001 \tapprox = \t{approxSecondDerivative(testFunction, pi / 4, 0.0000001)}")

h = 0.00001
xs = [0.8,0.9,1.0,1.1,1.2]
ys = []
for x in xs:
  ys.append(approxSecondDerivative(testFunction, x, h))
print(xs)
print(ys)
