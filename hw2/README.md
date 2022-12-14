# Tasksheet 2 by Kollin Murphy

## Task 1: Newton's Method

**Routine Name:** `newtons`

**Author:** Kollin Murphy

**Language:** Python3

**Description/Purpose**: The function `newtons` approximates the root of a function given the function, its first-order derivative, an initial approximation, the maximum number of iterations, and a tolerance.

#### Usage

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

#### Flags

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

#### Implementation

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

#### Driver Code

```python
from math import e

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


## Task 2: Secant Method

**Routine Name:** `secant`

**Author:** Kollin Murphy

**Language:** Python3

**Description/Purpose**: The function `secant` approximates the root of a function given the function, two initial approximations, the maximum number of iterations, and a tolerance level.

#### Usage

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

#### Flags

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

#### Implementation

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

#### Driver Code

```python
from math import e

def f(x):
  return x * pow(e, -x)

print(secant(f, 1, 0.5, 100, 0.00001))
```

Example output:

```
4.039302763706419e-07
```

## Task 3: Tabulation of Results

I modified the Newton's Method code as follows:

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

I modified the Secant code as follows:

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
    x2 = x0 - f0 / ((f1 - f0) / (x1 - x0))
    x0 = x1
    x1 = x2
    f0 = f1
    f1 = f(x1)
    i += 1

    if verbose:
      print("{:<12} {:<24} {:<24}".format(i, x1, abs(f1)))
  return x1
```



## Task 4: Hybrid Bisection / Newton's Method

**Routine Name:** `hybrid_newtons`

**Author**: Kollin Murphy

**Language**: Python3

**Description/Purpose**: The function `hybrid_newtons` approximates the root of a function given the function, its first-order derivative, initial lower and upper bounds, the maximum number of iterations, and a tolerance.

#### Usage

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

#### Flags

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

#### Implementation

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

#### Driver Code

```python
from math import cos, sin

def f(x):
  return 1000 * sin(x)

def fPrime(x):
  return 1000 * cos(x)

print(hybrid_newtons(f, fPrime, -2.1, 1, 100, 10e-12))
```

Example output:

```
Iteration    Algorithm    Approx. Root             Abs. Error              
0            bisection    -0.06562500000000003     0.19375                 
1            newtons      9.437033465323208e-05    0.06571937033465326     
2            newtons      -2.8014652501007903e-13  9.43703349333786e-05    
3            newtons      0.0                      2.8014652501007903e-13  
0.0
```


## Task 5: Hybrid Secant Method

**Routine Name:** `hybrid_secant`

**Author**: Kollin Murphy

**Language**: Python3

**Description/Purpose**: The function `hybrid_secant` approximates the root of a function given the function, two initial approximations, the maximum number of iterations, and a tolerance.

#### Usage

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

#### Flags

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

#### Implementation

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

#### Driver Code

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
Iteration    Algorithm    Approx. Root             Abs. Error              
0            bisection    -0.03125                 0.1875                  
1            secant       0.005670740864149072     0.03692074086414907     
2            secant       -7.555527984969349e-07   0.005671496416947569    
3            secant       4.048898900659039e-12    7.555568473958355e-07   
4            secant       -3.8531753143339277e-25  4.0488989006594245e-12  
-3.8531753143339277e-25
```
