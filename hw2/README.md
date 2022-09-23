# Tasksheet 2

See additional documentation [here](docs/README.md).

## Task 1: Newton's Method

```python
def newtons(f, fPrime, x, maxIterations, tolerance):
  i = 0
  y = 10.0 * tolerance
  while (abs(y) > tolerance and i < maxIterations):
    y = f(x) / fPrime(x)
    x = x - y
    i += 1
  return x
```

Driver code:

```python
from math import e

def f(x):
  return x * pow(e, -x)

def fPrime(x):
  return pow(e, -x) - (x * pow(e, -x))

print(newtons(f, fPrime, -3.0, 1000, 0.00001))
```



## Task 2: Secant Method

```python
def secant(f, x0, x1, maxIterations, tolerance):
  i = 0
  f0 = f(x0)
  f1 = f(x1)

  while (abs(f1) > tolerance and i < maxIterations):
    x2 = x0 - f0 / ((f1 - f0) / (x1 - x0))
    x0 = x1
    x1 = x2
    f0 = f1
    f1 = f(x1)
    i += 1
  return x1
```

Driver code:

```python
from math import e

def f(x):
  return x * pow(e, -x)

print(secant(f, 1, 0.5, 100, 0.00001))
```



## Task 3: Tabulation of Results

I modified the Newton's Method code as follows:

```python
from math import e
from sys import argv

verbose = False
for arg in argv:
  if arg == '-v' or arg == '--v' or arg == '--verbose':
    verbose = True

def newtons(f, fPrime, x, maxIterations, tolerance):
  i = 0
  y = 10.0 * tolerance

  if verbose:
    print("{:<12} {:<24} {:<24}".format("Iteration", "Approx. Root", "Abs. Error"))

  while (abs(y) > tolerance and i < maxIterations):
    y = f(x) / fPrime(x)
    x = x - y
    i += 1
    if verbose:
      print("{:<12} {:<24} {:<24}".format(i, x, abs(y)))
  return x

def f(x):
  return x * pow(e, -x)

def fPrime(x):
  return pow(e, -x) - (x * pow(e, -x))

print(newtons(f, fPrime, -3.0, 1000, 0.00001))
```

I modified the Secant code as follows:

```python
from math import e
from sys import argv

verbose = False
for arg in argv:
  if arg == '-v' or arg == '--v' or arg == '--verbose':
    verbose = True

def secant(f, x0, x1, maxIterations, tolerance):
  i = 0
  f0 = f(x0)
  f1 = f(x1)

  if verbose:
    print("{:<12} {:<24} {:<24}".format("Iteration", "Approx. Root", "Abs. Error"))

  while (abs(f1) > tolerance and i < maxIterations):
    x2 = x0 - f0 / ((f1 - f0) / (x1 - x0))
    x0 = x1
    x1 = x2
    f0 = f1
    f1 = f(x1)
    i += 1

    if verbose:
      print("{:<12} {:<24} {:<24}".format(i, x1, abs(f1)))
  return x1

def f(x):
  return x * pow(e, -x)

print(secant(f, 1, 0.5, 100, 0.00001))
```
