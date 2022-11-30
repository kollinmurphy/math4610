## Second Derivative Approximation

**Routine Name:** `approxSecondDerivative`

**Author**: Kollin Murphy

**Language**: Python3

**Description/Purpose**: The function `approxSecondDerivative` approximates the second derivative of a function given the function, a point at which to approximate the second derivative, and a step size.

#### Usage

The function accepts parameters as follows:

```python
approxSecondDerivative(function, x, stepSize)
```

The function returns an approximate second derivative of the given function.

Example usage:

```python
from math import pi, tan

def testFunction(x):
  return ((x - pi / 2) * tan(x) * tan(x)) / (x * x + 65)

print(f"h = 1 \t\tapprox = \t {approxSecondDerivative(testFunction, pi / 4, 1)}")
print(f"h = 0.1 \tapprox =\t{approxSecondDerivative(testFunction, pi / 4, 0.1)}")
print(f"h = 0.001 \tapprox = \t{approxSecondDerivative(testFunction, pi / 4, 0.001)}")
print(f"h = 0.00001 \tapprox = \t{approxSecondDerivative(testFunction, pi / 4, 0.00001)}")
print(f"h = 0.0000001 \tapprox = \t{approxSecondDerivative(testFunction, pi / 4, 0.0000001)}")
```

Example output:

```
h = 1           approx =         0.08888433309385914
h = 0.1         approx =        -0.0006870934574899935
h = 0.001       approx =        -6.767779469826762e-08
h = 0.00001     approx =        -6.767775576066448e-12
h = 0.0000001   approx =        -6.852157730108388e-16
```

#### Implementation

```python
def approxSecondDerivative(f, x0, h):
  return f(x0 - h) - 2 * f(x0) + f(x0 + h)
```