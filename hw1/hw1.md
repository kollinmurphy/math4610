# Math 4610 Task Sheet 1

**Kollin Murphy**

## Task 1

Solution (`rootFinder.py`):

```python
def fixedPointRootFinder(f, initial_approximation, tolerance, max_iterations):
  error = 10.0 * tolerance
  x0 = initial_approximation
  x1 = 0
  iterations = 0

  def g(x):
    return x - f(x)

  while error > tolerance and iterations < max_iterations:
    x1 = g(x0)
    error = abs(x1 - x0)
    x0 = x1
    iterations += 1
  return x0
```

Driver code (`driver.py`):

```python
from rootFinder import fixedPointRootFinder

def f(x):
  return x * pow(e, -x)

tolerance = 0.0000001
max_iterations = 1000
initial_guess = 0.5
print(fixedPointRootFinder(f, initial_guess, tolerance, max_iterations))

```

Example output:

```
4.676720855835697e-24
```

## Task 2

Solution (`rootFinder.py`):

```python
from math import e
from sys import argv

verbose = False
for arg in argv:
  if arg == '-v' or arg == '--v' or arg == '--verbose':
    verbose = True

def fixedPointRootFinder(f, initial_approximation, tolerance, max_iterations):
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

Driver code (`driver.py`):

```python
from rootFinder import fixedPointRootFinder

def f(x):
  return x * pow(e, -x)

tolerance = 0.0000001
max_iterations = 1000
initial_guess = 0.5
print(fixedPointRootFinder(f, initial_guess, tolerance, max_iterations))
```

To run in terminal:

```bash
% python driver.py -v
```

Example output:

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



## Task 3

Given function:

$$
\begin{align*}
	f(x) = 10.14 * e^{x^2} \cos(\frac{\pi}{x})
\end{align*}
$$

I wrote a program to check every 0.25 along the interval [-3, 7].

```python
from math import cos, e, pi
from rootFinder import fixedPointRootFinder

def f(x):
  return 10.14 * pow(e, x * x) * cos(pi / x)

for i in range(-12, 29):
  print()
  print(i/4)
  try:
    print(fixedPointRootFinder(f, i/4, 0.01, 100))
  except:
    print('failure')
```

I tested this using `g(x) = x - f(x)` and also attempting it with `g(x) = x + f(x)`. Every point failed except at `x=-2` and `x=2`, which produced approximate results of `-2.0000000000000338` and `1.999999999999966`, respectively.

## Task 4

Solution:

```python
from math import ceil, e, log
from sys import argv

verbose = False
for arg in argv:
  if arg == '-v' or arg == '--v' or arg == '--verbose':
    verbose = True

def bisect(f, a, b, tol):
  fa = f(a)
  fb = f(b)

  if fa * fb >= 0:
    raise Exception('Intermediate Value Theorem is not satisfied by the initial conditions')

  k = ceil(-log(tol / (b - a)) / log(2))
  c = a

  if verbose:
    print("{:<12} {:<24} {:<24}".format("Iteration", "Approx. Root", "Abs. Error"))

  for i in range(0, k):
    c = (a + b) / 2
    fc = f(c)
    if fa * fc < 0:
      b = c
      fb = fc
    else:
      a = c
      fa = fc
    if verbose:
      error = abs(a - b)
      print("{:<12} {:<24} {:<24}".format(i, c, error))
  
  return c
```

###### Test case (`f` from Task 1)

Driver code:

```python
def f(x):
  return x * pow(e, -x)

print(bisect(f, -1, 2, 0.0001))
```

Output:

```
Iteration    Approx. Root             Abs. Error              
0            0.5                      1.5                     
1            -0.25                    0.75                    
2            0.125                    0.375                   
3            -0.0625                  0.1875                  
4            0.03125                  0.09375                 
5            -0.015625                0.046875                
6            0.0078125                0.0234375               
7            -0.00390625              0.01171875              
8            0.001953125              0.005859375             
9            -0.0009765625            0.0029296875            
10           0.00048828125            0.00146484375           
11           -0.000244140625          0.000732421875          
12           0.0001220703125          0.0003662109375         
13           -6.103515625e-05         0.00018310546875        
14           3.0517578125e-05         9.1552734375e-05        
3.0517578125e-05
```

###### Test case (`f` from Task 3)

Driver code:

```python
def f(x):
  return 10.14 * pow(e, x * x) * cos(pi / x)

print(bisect(f, 1.8, 2.1, 0.0001))
```

Output:

```
Iteration    Approx. Root             Abs. Error              
0            1.9500000000000002       0.1499999999999999      
1            2.0250000000000004       0.07500000000000018     
2            1.9875000000000003       0.03750000000000009     
3            2.0062500000000005       0.018750000000000266    
4            1.9968750000000004       0.009375000000000133    
5            2.0015625000000004       0.004687499999999956    
6            1.9992187500000003       0.002343750000000089    
7            2.0003906250000005       0.0011718750000002665   
8            1.9998046875000004       0.0005859375000001332   
9            2.0000976562500004       0.0002929687499999556   
10           1.9999511718750003       0.00014648437500008882  
11           2.0000244140625005       7.324218750026645e-05   
2.0000244140625005
```

## Task 5

You can find my GitHub repository for this course [here](https://github.com/kollinmurphy/math4610).
