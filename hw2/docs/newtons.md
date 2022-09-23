# Newton's Method

**Routine Name:** `newtons`

**Author:** Kollin Murphy

**Language:** Python3

**Description/Purpose**: The function `newtons` approximates the root of a function given the function, its first-order derivative, an initial approximation, the maximum number of iterations, and a tolerance.

## Usage

The function accepts parameters as follows:

```python
newtons(function, functionDerivative, initialApproximation, maximumIterations, tolerance)
```

Example usage:

```python
def f(x):
  return x * pow(e, -x)

def fPrime(x):
  return pow(e, -x) - (x * pow(e, -x))

print(newtons(f, fPrime, -3.0, 1000, 0.00001))
```

Example output:

```
-1.3178467639212633e-14
```

## Flags

This command also accepts a verbose flag that is given from the command line with `-v`, `--v`, or `--verbose`.

When verbose mode is enabled, it will output a table of results. From the same example above, it would produce the following table.

```
Iteration    Approx. Root             Abs. Error              
1            -2.25                    0.75                    
2            -1.5576923076923077      0.6923076923076923      
3            -0.9486697513013304      0.6090225563909774      
4            -0.461840338227204       0.48682941307412636     
5            -0.14590957195262932     0.31593076627457467     
6            -0.018578781178275744    0.12733079077435358     
7            -0.00033887522148356997  0.018239905956792174    
8            -1.1479751370265271e-07  0.0003387604239698673   
9            -1.3178467639212633e-14  1.1479750052418507e-07  
```

## Implementation

`newtons` is implemented as follows:

```python
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
```
