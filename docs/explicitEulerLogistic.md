
## Explicit Euler Logistic

**Routine Name:** `explicitEulerLogistic`

**Author**: Kollin Murphy

**Language**: Python3

**Description/Purpose**: The function `explicitEulerLogistic` approximates the solution to the logistic differential equation using the explicit Euler method.

#### Usage

The function accepts parameters as follows:

```python
explicitEulerLogistic(alpha, beta, p0)
```

It returns a tuple of x values and y values.

Example usage:

```python
print(explicitEulerLogistic(0.2, 0.0005, 10.0))
```

Example output:

```
([0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0], [10, 11.95, 14.268598749999999, 17.02052204485575, 20.27977736848719, 24.130098157126923, 28.664986970016017, 33.98714362502363, 40.20700938413436, 47.4401094591533, 55.80284935823573, 65.40644023163391, 76.34872706607354, 88.70390841698335, 102.51049841615574, 117.75839695662256, 134.37655632105032, 152.22333814090817, 171.08203343170916, 190.66390903648676, 210.62032773924727, 230.56393205860275, 250.09685508716132, 268.8420076423493, 286.47239663423466, 302.7336589444005, 317.4565566043483, 330.5585352596731, 342.03576969509743, 351.9486897586581, 360.4044875989716, 367.5396877780273, 373.50491428764786, 378.4529366466659, 382.5302113477564, 385.8715723204281, 388.597451621994, 390.81295224283883, 392.6081608710249, 394.05920905396556, 395.2297207446402, 396.1723988136251, 396.9305937854791, 397.5397644011284, 398.0287851413018, 398.4210852690319, 398.7356217293617, 398.98769805728375, 399.1896460682154, 399.35138851782494, 399.4809004658325, 399.58458564050284, 399.6675822278572, 399.73401053149814, 399.78717304999986, 399.8297157923446, 399.86375813551996, 399.89099722749313, 399.9127918411923, 399.9302296703224, 399.94418130230844, 399.95534348398326, 399.9642737900844, 399.9714183938865, 399.9771343066551, 399.9817071839041, 399.9853655798097, 399.98829235676465, 399.9906338168773, 399.9925070096391, 399.99400557963884, 399.99520444574455, 399.996163545097, 399.9969308287184, 399.9975446582648, 399.99803572359747, 399.9984285769488, 399.9987428603244, 399.9989942874693, 399.99919542946975, 399.99935634325215, 399.99948507439456, 399.9995880593831, 399.99967044742164, 399.99973635788297, 399.9997890862716, 399.9998312689951, 399.99986501518185, 399.99989201213634, 399.99991360970324, 399.99993088775886, 399.9999447102047, 399.99995576816224, 399.99996461452884, 399.9999716916225, 399.9999773532976, 399.9999818826378, 399.9999855061101, 399.999988404888, 399.9999907239103])
```

Plot of reults:

![image-20221016161651267](../hw4/img/task5-test1.png)

#### Implementation

```python
def explicitEulerLogistic(alpha, beta, p0):
  def f(p):
    return alpha * p - beta * (p * p)

  T = 100
  n = 100
  t0 = 0
  h = (T - t0) / n
  f0 = f(p0)
  tvals = [t0]
  pvals = [p0]

  for i in range(1, n):
    t1 = t0 + h
    p1 = p0 + h * f0
    
    f0 = f(p1)
    t0 = t1
    p0 = p1

    tvals.append(t1)
    pvals.append(p1)
  return (tvals, pvals)
```