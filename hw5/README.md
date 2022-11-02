# Tasksheet 5 by Kollin Murphy

## Task 1

I implemented the Implicit Euler Method to solve the logistics equation. I used the following code to solve the logistics equation:

```python
def implicitEulerMethod(alpha, beta, p0, n=100, h=1.0):
  """Solve the logistic equation using the implicit Euler method."""

  def fPrime(p):
    return 1 - alpha * h + 2 * h * beta * p

  ptable = [0 for i in range(n)]
  ptable[0] = p0
  for i in range(1,n):
    def f(p):
      return p - ptable[i - 1] - h * (alpha * p - beta * p * p)
    pval = newtons(f, fPrime, ptable[i - 1], 100, 1e-6)
    ptable[i] = pval
  return zip(*[[i*h, ptable[i]] for i in range(n)])
```

This function takes as input the parameters alpha, beta, p0, n, and h. It returns two lists, the first one containing the values of $t$, and the second one containing the values of $P(t)$. The function uses the implicit Euler method to solve the logistics equation. The function uses Newton's method to solve the rootfinding problem that arises during the implicit Euler algorithm. 

I used the following code to test and plot the results, and save the plots as individual files.

```python
tests = [
  [0.2, 0.0005, 10.0, 75, 1.0, "Test 1", "task1-test1.png"],
  [0.01, 0.0005, 10.0, 350, 2.0, "Test 2", "task1-test2.png"],
  [2.0, 0.0005, 10.0, 75, 0.1, "Test 3", "task1-test3.png"],
]

for test in tests:
  tvals, pvals = implicitEulerMethod(test[0], test[1], test[2], test[3], test[4])
  plt.title(test[5])
  plt.xlabel('t')
  plt.ylabel('P(t)')
  plt.plot(tvals, pvals)
  plt.savefig(test[6])
  plt.clf()
```

I varied the values of $n$ and $h$ to fit the curves in a window where we could see them reach their carrying capacity. These are my results:

**Test 1.** I used the following parameters: $\alpha = 0.2$, $\beta = 0.0005$, $P_0 = 10.0$, $n = 75$, $h = 1.0$. This is the plot:

<img src="task1-test1.png" width="400px" />

**Test 2.** I used the following parameters: $\alpha = 0.01$, $\beta = 0.0005$, $P_0 = 10.0$, $n = 350$, $h = 2.0$. This is the plot:

<img src="task1-test2.png" width="400px" />

**Test 3.** I used the following parameters: $\alpha = 2.0$, $\beta = 0.0005$, $P_0 = 10.0$, $n = 75$, $h = 0.1$. This is the plot:

<img src="task1-test3.png" width="400px" />

## Task 2

The given differential equation is:
$$
\begin{align*}
\frac{dP}{dt} &= \alpha P - \beta P^2
\end{align*}
$$
Below is a complete analytic solution to the differential equation.
$$
\begin{align*}
\frac{dP}{dt} &= \alpha P - \beta P^2 \\
dP &= (\alpha P - \beta P^2)\;dt \\
dt &= \frac{dP}{\alpha P - \beta P^2} \\
\int 1 \;dt &=\int \frac{1}{\alpha  P - \beta P^2}\;dP \\
t+C &= \int \frac{1}{P^2(\frac {\alpha}{P}  - \beta)}\;dP
\end{align*}
$$
By using $u$-substitution, we are able to solve this integral.
$$
\begin{align*}
u &= \frac{\alpha}{P} - \beta \\
du &=-\alpha P^{-2}\;dP
\end{align*}
$$
Substituting $u$ and $du$ into the integral, we get:
$$
\begin{align*}
\int \frac{1}{P^2(\frac {\alpha}{P}  - \beta)}\;dP &= \int -\frac{1}{\alpha u}\;du \\
&=-\frac{1}{\alpha}\ln|u|+C \\
&= -\frac{\ln|\frac{\alpha}{P} - \beta |}{\alpha} + C
\end{align*}
$$
Therefore,
$$
\begin{align*}
t &= -\frac{\ln|\frac{\alpha}{P} - \beta |}{\alpha} + C
\end{align*}
$$

Then, we solve for $t$:
$$
\begin{align*}
-\alpha t+c_1 &= \ln |\frac{\alpha}{P} - \beta | \\
e^{-\alpha t+c_1} &=\frac{\alpha}{P} - \beta \\
c_2e^{-\alpha t}+\beta &=\frac{\alpha}{P} \\
P &= \frac{\alpha}{c_2e^{-\alpha t}+\beta}
\end{align*}
$$
We can also solve for $c_2$ by knowing that $P(0)=P_0$.
$$
\begin{align*}
P(0) = P_0 &= \frac{\alpha}{c_2e^{-\alpha (0)}+\beta} \\
P_0 &= \frac{\alpha}{c_2 + \beta} \\
c_2 &= \frac{\alpha}{P_0}-\beta
\end{align*}
$$
So, the final analytic solution is $P(t, \alpha, \beta, P_0) = \frac{\alpha}{(\frac{\alpha}{P_0}-\beta)e^{-\alpha t}+\beta}$.

I was able to implement this analytic solution in Python using the following code:

```python
from math import exp

def analyticSolution(alpha, beta, p0, n=100, h=1.0):
  """Solve the logistic equation using the analytic solution."""
  c = alpha / p0 - beta
  return zip(*[[i*h, alpha / (c * exp(-alpha * i*h) + beta)] for i in range(n)])
```

This method returns two lists, the first one containing the values of $t$, and the second one containing the values of $P(t)$.

## Task 3

In this task, I used the same test cases as in Task 1. I used the following code to test and plot the results, and save the plots as individual files.

```python
tests = [
  [0.2, 0.0005, 10.0, 50, 2, "Test 1 (n=50, h=2.0)", "task3-test1-a.png"],
  [0.2, 0.0005, 10.0, 200, 0.5, "Test 1 (n=200, h=0.5)", "task3-test1-b.png"],
  [0.01, 0.0005, 10.0, 50, 15.0, "Test 2 (n=50, h=15.0)", "task3-test2-a.png"],
  [0.01, 0.0005, 10.0, 200, 3.75, "Test 2 (n=200, h=3.75)", "task3-test2-b.png"],
  [2.0, 0.0005, 10.0, 50, 0.2, "Test 3 (n=50, n=0.2)", "task3-test3-a.png"],
  [2.0, 0.0005, 10.0, 200, 0.05, "Test 3 (n=200, 0.05)", "task3-test3-b.png"],
]

for test in tests:
  tImplicit, pImplicit = implicitEulerMethod(test[0], test[1], test[2], test[3], test[4])
  tExplicit, pExplicit = explicitEulerLogistic(test[0], test[1], test[2], test[3], test[4])
  tAnalytic, pAanalytic = analyticSolution(test[0], test[1], test[2], test[3], test[4])

  plt.title(test[5])
  plt.xlabel('t')
  plt.ylabel('P(t)')

  plt.plot(tImplicit, pImplicit, 'r', label='Implicit Euler')
  plt.plot(tExplicit, pExplicit, 'g', label='Explicit Euler')
  plt.plot(tAnalytic, pAanalytic, 'b', label='Analytic')

  plt.legend(loc="lower right")
  plt.savefig(test[6])
  plt.clf()
```

For each test, I plotted the results of the implicit Euler method, the explicit Euler method, and the analytic solution. I also saved the plots as individual files. I used two different values for $n$ for each test case, $n=50$ and $n=200$. I varied the value of $h$ to fit the plot in a window that would show the entire relevant portion of the graph, and to display the same view of the graph with the different values of $n$. Below are the plots for each test case.

**Test 1.** I used the following parameters: $\alpha = 0.2$, $\beta = 0.0005$, $P_0 = 10.0$.

With $n=50$ and $h=2.0$:

<img src="task3-test1-a.png" width="400px" />

With $n=200$ and $h=0.5$:

<img src="task3-test1-b.png" width="400px" />

**Test 2.** I used the following parameters: $\alpha = 0.01$, $\beta = 0.0005$, $P_0 = 10.0$. This is the plot:

With $n=50$ and $h=15.0$:

<img src="task3-test2-a.png" width="400px" />

With $n=200$ and $h=3.75$:

<img src="task3-test2-b.png" width="400px" />

**Test 3.** I used the following parameters: $\alpha = 2.0$, $\beta = 0.0005$, $P_0 = 10.0$. This is the plot:

With $n=50$ and $h=0.2$:

<img src="task3-test3-a.png" width="400px" />

With $n=200$ and $h=0.05$:

<img src="task3-test3-b.png" width="400px" />

As seen from the graphs, there is a little deviation from the analytic solution for both the implicit and explicit methods of approximation. This is expected, since the analytic solution is exact, while the numerical methods are approximations. The deviation was larger when the value of $h$ was larger, and smaller when the value of $h$ was smaller. This is also expected, since the smaller the value of $h$, the more accurate the approximation.

## Task 4

I implemented the trapezoid rule using Python. I used the following code to implement the trapezoid rule:

```python
def trapezoidRule(f, a, b, n):
  h = (b - a) / n
  sum = 0.5 * (f(a) + f(b))
  for i in range(1, n):
    sum += f(a + i * h)
  return sum * h
```

This function takes as input a function $f$, the interval $[a, b]$, and the number of subintervals $n$. It returns the approximation of the integral of $f$ over the interval $[a, b]$ using the composite trapezoid rule.

I tested the function by using the following code:

```python
from math import e, pow, pi

tests = [2, 4, 8, 16, 100000]

def f(x):
  return pow(e, -x * x)

for n in tests:
  print(f"n={n}; approx={trapezoidRule(f, 0, pi / 4, n)}")
```

The output of the test code is:

```
n=2; approx=0.6388862805734845
n=4; approx=0.6471507696813964
n=8; approx=0.6491991053630145
n=16; approx=0.6497100964398593
n=1000; approx=0.6498802865050275
```

It appears that the sequence is converging to a value of approximately `0.64988`. This is the correct value of the integral of `e^(-x^2)` from `0` to `pi / 4`.

## Task 5

I implemented the Simpson's rule using Python. I used the following code to implement the Simpson's rule:

```python
def simpsonsRule(f, a, b, n):
  h = (b - a) / n
  x = a
  sum = 0
  for i in range(n):
    sum += f(x) + 4 * f(x + h / 2) + f(x + h)
    x += h
  return sum * h / 6
```

The function takes as input a function $f$, the interval $[a, b]$, and the number of subintervals $n$. It returns the approximation of the integral of $f$ over the interval $[a, b]$ using the composite Simpson's rule.

I tested my implementation with the following code:

```python
from math import e, pow, pi

tests = [2, 4, 8, 16, 100000]

def f(x):
  return pow(e, -x * x)

for test in tests:
  print(f"n={test}; approx={simpsonsRule(f, 0, pi / 4, test)}")
```

The output of the test code is:

```
n=2; approx=0.6499055993840338
n=4; approx=0.6498818839235537
n=8; approx=0.6498804267988076
n=16; approx=0.6498803361175072
n=100000; approx=0.6498803300785783
```

It appears that the sequence is converging to a value of approximately `0.64988`. This is the correct value of the integral of `e^(-x^2)` from `0` to `pi / 4`. Using Wolfram Alpha, I found that the exact value of the integral is `0.6498803300786573037276521829129935239240253152760926685227030955`. This is very close to the value I found using Simpson's rule. I used this value to perform a convergence analysis on my code. I used the following code to perform the convergence analysis:

```python
def calculateNFromH(a, b, h):
  return (b - a) / h

hvals = [1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4196]
nvals = [int(calculateNFromH(0, 1, h)) for h in hvals]

# I calculated the actual value of the integral using Wolfram Alpha
actual = 0.6498803300786573037276521829129935239240253152760926685227030955

def f(x):
  return pow(e, -x * x)

approximations = [simpsonsRule(f, 0, pi / 4, n) for n in nvals]
errors = [abs(actual - approximation) for approximation in approximations]

plt.plot([log(h) for h in hvals], [log(error) for error in errors])
plt.title("Log-Log Plot of Error vs. h")
plt.xlabel("log(h)")
plt.ylabel("log(error)")
plt.savefig("task5.png")
plt.show()
```

The resulting plot is:

<img src="task5.png" width="400px" />

As can be seen in the plot, as the value of h decreased, the accuracy increased until a certain point at which the accuracy stopped increasing and actually decreased due to floating point errors.

The values returned from the test code are:

```
h              | approx         | error
------------------------------------------------
2.50000000e-01 | 6.49881884e-01 | 1.55384490e-06
1.25000000e-01 | 6.49880427e-01 | 9.67201503e-08
6.25000000e-02 | 6.49880336e-01 | 6.03884986e-09
3.12500000e-02 | 6.49880330e-01 | 3.77331721e-10
1.56250000e-02 | 6.49880330e-01 | 2.35813591e-11
7.81250000e-03 | 6.49880330e-01 | 1.47415413e-12
3.90625000e-03 | 6.49880330e-01 | 9.24815780e-14
1.95312500e-03 | 6.49880330e-01 | 5.10702591e-15
9.76562500e-04 | 6.49880330e-01 | 9.99200722e-16
4.88281250e-04 | 6.49880330e-01 | 8.88178420e-16
2.38322212e-04 | 6.49880330e-01 | 9.99200722e-15
```

When $h$ was cut in half, the accuracy of the approximation increased by a factor of 16. This is because the number of intervals was doubled. This is consistent with the convergence rate of Simpson's rule, which is $O(h^4)$. This means that the error is proportional to $h^4$.
