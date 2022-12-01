# Tasksheet 9 by Kollin Murphy

## Task 1: Inverse Power Method

I implemented the inverse power method to find the smallest eigenvalue of a matrix. I used Jacobi iteration to solve the system of linear equations that arose during the algorithm. My implementation is as follows.

```c
double inversePowerMethod(int n, double *A, double *result, double *x0, double tolerance, int maxIterations)
{
  double lambda0 = 0.0;
  double lambda1 = 0.0;
  double error = 10.0 * tolerance;
  int iter = 0;
  double y[n];
  double z[n];
  double w[n];
  double initialForJacobi[n];
  double norm;

  for (int i = 0; i < n; i++)
  {
    initialForJacobi[i] = 1.0;
    y[i] = 0.0;
    w[i] = 0.0;
  }

  jacobiIteration(n, A, x0, initialForJacobi, y, tolerance, maxIterations);

  while (error > tolerance && iter < maxIterations)
  {
    norm = l2Norm(n, y);
    for (int i = 0; i < n; i++)
      z[i] = y[i] / norm;

    jacobiIteration(n, A, z, initialForJacobi, w, tolerance, maxIterations);
    lambda1 = dotProduct(n, z, w);
    error = fabs(lambda1 - lambda0);
    lambda0 = lambda1;
    iter++;
    for (int i = 0; i < n; i++)
      y[i] = w[i];
  }

  for (int i = 0; i < n; i++)
    result[i] = y[i];
  return 1.0 / lambda1;
}
```


To test the algorithm, I used the matrix:

$$
A = \begin{bmatrix}
67 & 8 & 4 \\
8 & 90 & 8 \\
4 & 8 & 83
\end{bmatrix}
$$

I used the parameters $x_0 = (1, 1, 1)$, $tolerance = 10^{-8}$, and $maxIterations = 100$. The result was:

```
lambda = 64.363547
v = [0.014888, -0.004228, -0.001363]
```

This is the correct result, as the smallest eigenvalue of $A$ is $64.363547$ and the corresponding eigenvector is $[0.014888, -0.004228, -0.001363]$.

## Task 2: Shifted Inverse Power Method

