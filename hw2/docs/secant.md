# Secant Method

**Routine Name:** `secant`

**Author:** Kollin Murphy

**Language:** Python3

**Description/Purpose**: The function `secant` approximates the root of a function given the function, two initial approximations, the maximum number of iterations, and a tolerance level.

## Usage

The function accepts parameters as follows, where `x0` and `x1` are two initial approximations.

```python
secant(function, x0, x1, maximumIterations, tolerance)
```

Example usage:

```python
def f(x):
  return x * pow(e, -x)

print(secant(f, 1, 0.5, 100, 0.00001))
```

Example output:

```
4.039302763706419e-07
```

## Flags

This command also accepts a verbose flag that is given from the command line with `-v`, `--v`, or `--verbose`.

When verbose mode is enabled, it will output a table of results. From the same example above, it would produce the following table.

```
Iteration    Approx. Root             Abs. Error              
1            -1.8467422493615944      11.706747557360922
2            0.44074231484814863      0.2836434721461853
3            0.3866298123720783       0.2626540759096701
4            -0.2905154076645947      0.3884529365322745
5            0.11347202321181965      0.10129978529945151
6            0.029911813768984136     0.029030346036369014
7            -0.003653989698321175    0.00366736576215364
8            0.00011074097071405798   0.00011072870783047764
9            4.039302763705064e-07    4.039301132108712e-07
```

## Implementation

`secant` is implemented as follows:

```python
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
    x2 = x1 - f1 * (x1 - x0)  /(f1 - f0)
    x0 = x1
    x1 = x2
    f0 = f1
    f1 = f(x1)
    i += 1

    if verbose:
      print("{:<12} {:<24} {:<24}".format(i, x1, abs(f1)))
  return x1
```
