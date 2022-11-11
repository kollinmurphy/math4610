# Tasksheet 6 by Kollin Murphy

## Task 1: Computing an Approximation of Pi

To calculate an approximation of $\pi$, I decided to use the following integral:
$$
\int \frac{1}{\sqrt{1-x^2}}dx = \arcsin(x)+C
$$
 By taking the definite integral from $0$ to $1$, we are able to get a multiple of $\pi$.
$$
\begin{align*}
\int_0^{\frac{1}{2}} \frac{1}{\sqrt{1-x^2}}dx &= \arcsin(\frac{1}{2})-\arcsin(0) \\
&= \frac{\pi}{6}-0 \\
&= \frac{\pi}{6}
\end{align*}
$$
Therefore, we can use the following to generate our approximation:
$$
6\int_0^{\frac{1}{2}} \frac{1}{\sqrt{1-x^2}}dx=\pi
$$

## Task 2: Coding an Approximation of Pi

To code an approximation of $\pi$, I used the trapezoid rule to approximate the integral. I used the following code to do this:

```c
#include <stdio.h>
#include <math.h>

double f(double x) {
  return 1 / sqrtf(1 - (x * x));
}

double trapezoidRule(double (*f)(double), double a, double b, int n) {
  double h = (b - a) / n;
  double sum = 0.5 * (f(a) + f(b));
  for (int i = 1; i < n; ++i) {
    sum += f(a + i * h);
  }
  return sum * h;
}

int main() {
  double pi = 6 * trapezoidRule(f, 0, 0.5, 100000000);
  printf("pi = %.16f\n", pi);
  return 0;
}
```

This code is serial. It uses a total of $100000000$ iterations to approximate $\pi$. The output of this code is:

```
pi = 3.1415926573648885
```

It ran in $0.201$ seconds on my machine as timed by using the `time` command.

## Task 3: Parallelizing an Approximation of Pi

I added OpenMP directives to parallelize the code. This is my implementation:

```c
#include <stdio.h>
#include <math.h>
#include <omp.h>

double f(double x) {
  return 1 / sqrtf(1 - (x * x));
}

static long num_steps = 100000000;
#define NUM_THREADS 4

double trapezoidRule(double (*f)(double), double a, double b, int n) {
  double h = (b - a) / n;
  double sum[NUM_THREADS];
  int actualNumThreads = 0;
  double accumulator = 0.5 * (f(a) + f(b));

  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numThreads = omp_get_num_threads();
    if (id == 0) actualNumThreads = numThreads;

    for (int i = 1 + id; i < n; i = i + numThreads) {
      sum[id] += f(a + i * h);
    }
  }

  for (int i = 0; i < actualNumThreads; ++i) {
    accumulator += sum[i];
  }
  return accumulator * h;
}

int main() {
  double pi = 6 * trapezoidRule(f, 0, 0.5, num_steps);
  printf("pi = %.16f\n", pi);
  return 0;
}
```

As before, I got the same output `pi = 3.1415926573648885`. With the compiler optimizations and set to run on 4 threads, it ran in $0.174$ seconds, a little bit faster than the serial code but definitely not as much of a decrease as I expected. I assume that this is due to increased overhead by setting up and taking down the threads. The calculations done here just don't take long enough to start seeing massive benefits from parallelization.

## Task 5: Linear Algebra Algorithms

I implemented each algorithm. Below is the implementation of each function.

```python
TODO
```

The following is the test function I used to test the algorithms.

```python
TODO
```

The following is the output of the test function.

```
TODO
```
