# Hybrid Bisection / Secant Method

**Routine Name:** `hybrid_secant`

**Author**: Kollin Murphy

**Language**: Python3

**Description/Purpose**: The function `hybrid_secant` approximates the root of a function given the function, two initial approximations, the maximum number of iterations, and a tolerance.

## Usage

The function accepts parameters as follows:

```python
hybrid_secant(function, x0, x1, maximumIterations, tolerance)
```

Example usage:

```python
from math import cos, sin

def f(x):
    return 1000 * sin(x)

def fPrime(x):
    return 1000 * cos(x)

print(hybrid_secant(f, -2, 1, 100, 10e-12))
```

Example output:

```
-3.8531753143339277e-25
```

## Flags

This command also accepts a verbose flag that is given from the command line with `-v`, `--v`, or `--verbose`.

When verbose mode is enabled, it will output a table of results. A sample output table may look similar to the following.

```
Iteration    Algorithm    Approx. Root             Abs. Error              
0            bisection    -0.03125                 0.1875                  
1            secant       0.005670740864149072     0.03692074086414907     
2            secant       -7.555527984969349e-07   0.005671496416947569    
3            secant       4.048898900659039e-12    7.555568473958355e-07   
4            secant       -3.8531753143339277e-25  4.0488989006594245e-12
```

## Implementation

`hybrid_secant` is implemented as follows:

```python
from sys import argv

verbose = False
for arg in argv:
    if arg == '-v' or arg == '--v' or arg == '--verbose':
        verbose = True

bisection_iterations = 4

def hybrid_secant(f, x0, x1, maxIterations, tolerance):
    lastA = x0
    lastB = x1
    error = 10.0 * tolerance
    i = 0

    if verbose:
        print("{:<12} {:<12} {:<24} {:<24}".format("Iteration", "Algorithm", "Approx. Root", "Abs. Error"))

    while (error > tolerance and i < maxIterations):
        f0 = f(x0)
        f1 = f(x1)
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        e = abs(x2 - x1)

        if (e > error):
            a = x0
            b = x1
            fa = f(a)
            fb = f(b)
            if fa * fb > 0:
                a = lastA
                b = lastB
                fa = f(a)
                fb = f(b)
                if fa * fb > 0:
                    raise Exception('Intermediate Value Theorem is not satisfied by the bisection endpoints')
            for _ in range(bisection_iterations):
                c = (a + b) / 2
                fc = f(c)
                if fa * fc < 0:
                    b = c
                    fb = fc
                else:
                    a = c
                    fa = fc
            x0 = x1
            x1 = (a + b) / 2
            error = abs(a - b)
            lastA = a
            lastB = b
            if verbose:
                print("{:<12} {:<12} {:<24} {:<24}".format(
                    i, 'bisection', x1, error))
        else:
            x0 = x1
            x1 = x2
            error = e
            if verbose:
                print("{:<12} {:<12} {:<24} {:<24}".format(
                    i, 'secant', x1, error))
        i += 1
    return x1
```

