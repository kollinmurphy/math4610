## Fixed Point Iteration

**Routine Name:** `fixedPoint`

**Author:** Kollin Murphy

**Language:** Python3

**Description/Purpose**: The function `fixedPoint` approximates the root of a function given the function, an initial approximation, a tolerance, and the maximum number of iterations. It uses the fixed point iteration algorithm to do so.

#### Usage

The function accepts parameters as follows:

```python
fixedPoint(f, initial_approximation, tolerance, max_iterations)
```

The function returns the approximate root of the function `f` as a float.

Example usage:

```python
from rootFinding import fixedPoint

def f(x):
  return x * pow(e, -x)

print(fixedPoint(f, initial_approximation=0.5, tolerance=0.0000001, max_iterations=1000))
```

Example output:

```
4.676720855835697e-24
```

#### Flags

This command also accepts a verbose flag that is given from the command line with `-v`, `--v`, or `--verbose`.

When verbose mode is enabled, it will output a table of results. From the same example above, it would produce the following table.

```
Iteration    Approx. Root             Abs. Error              
1            0.1967346701436833       0.3032653298563167      
2            0.035135130278937304     0.16159953986474598     
3            0.0012130423916044478    0.03392208788733286     
4            1.4705797257291708e-06   0.0012115718118787186   
5            2.162603139558558e-12    1.4705775631260313e-06  
6            4.676720855835697e-24    2.1626031395538813e-12  
4.676720855835697e-24
```

#### Implementation

`fixedPoint` is implemented as follows:

```python
from math import e
from sys import argv

verbose = False
for arg in argv:
  if arg == '-v' or arg == '--v' or arg == '--verbose':
    verbose = True

def fixedPoint(f, initial_approximation, tolerance, max_iterations):
  error = 10.0 * tolerance
  x0 = initial_approximation
  x1 = 0
  iterations = 0

  def g(x):
    return x - f(x)
  
  if verbose:
    print("{:<12} {:<24} {:<24}".format("Iteration", "Approx. Root", "Abs. Error"))

  while error > tolerance and iterations < max_iterations:
    x1 = g(x0)
    error = abs(x1 - x0)
    x0 = x1
    iterations += 1

    if verbose:
      print("{:<12} {:<24} {:<24}".format(iterations, x0, error))
  return x0
```

