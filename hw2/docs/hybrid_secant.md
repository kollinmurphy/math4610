# Hybrid Bisection / Secant Method

**Routine Name:** `hybrid_newtons`

**Author**: Kollin Murphy

**Language**: Python3

**Description/Purpose**: The function `hybrid_newtons` approximates the root of a function given the function, its first-order derivative, initial lower and upper bounds, the maximum number of iterations, and a tolerance.

## Usage

The function accepts parameters as follows:

```python
hybrid_newtons(function, functionDerivative, lower, upper, maximumIterations, tolerance)
```

Example usage:

```python
from math import cos, e, pi, sin

def f(x):
  return 10.14 * pow(e, pow(x, 2)) * cos(pi / x)

def fPrime(x):
  return ((1014 * cos(pi / x) * pow(x, 3)) + (507 * pi * sin(pi / x))) * pow(e, pow(x, 2)) / (50 * pow(x, 2))

print(hybrid_newtons(f, fPrime, -3, 7, 10, 0.001))
```

Example output:

```
2.0
```

## Flags

This command also accepts a verbose flag that is given from the command line with `-v`, `--v`, or `--verbose`.

When verbose mode is enabled, it will output a table of results. A sample output table may look similar to the following.

```
Iteration    Algorithm    Approx. Root             Abs. Error              
0            bisection    -0.06562500000000003     0.19375                 
1            newtons      9.437033465323208e-05    0.06571937033465326     
2            newtons      -2.8014652501007903e-13  9.43703349333786e-05    
3            newtons      0.0                      2.8014652501007903e-13  
0.0
```

## Implementation

`hybrid_newtons` is implemented as follows:

```python
from sys import argv

verbose = False
for arg in argv:
  if arg == '-v' or arg == '--v' or arg == '--verbose':
    verbose = True

bisection_iterations = 4

def hybrid_newtons(f, fPrime, a, b, maxIterations, tolerance):
  error = 10.0 * tolerance
  i = 0

  if verbose:
    print("{:<12} {:<12} {:<24} {:<24}".format("Iteration", "Algorithm", "Approx. Root", "Abs. Error"))

  x0 = (a + b) / 2
  while (error > tolerance and i < maxIterations):
    x1 = x0 - (f(x0) / fPrime(x0))
    e1 = abs(x1 - x0)

    if (e1 > error):
      fa = f(a)
      for _ in range(bisection_iterations):
        c = (a + b) / 2
        fc = f(c)
        if fa * fc < 0:
          b = c
          fb = fc
        else:
          a = c
          fa = fc
      x0 = (a + b) / 2
      error = abs(a - b)
      if verbose:
        print("{:<12} {:<12} {:<24} {:<24}".format(i, 'bisection', x0, error))
    else:
      x0 = x1
      error = e1
      if verbose:
        print("{:<12} {:<12} {:<24} {:<24}".format(i, 'newtons', x0, error))
    i += 1
  return x0
```

