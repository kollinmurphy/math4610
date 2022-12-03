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

I used the same algorithm as in Task 1, the Inverse Power Method. The only difference was that I subtracted $\mu$ from the diagonal of the matrix before running the algorithm. I used the same matrix as in Task 1, but with $\mu = 77$. The result was:

$$
A = \begin{bmatrix}
-10 & 8 & 4 \\
8 & 13 & 8 \\
4 & 8 & 6
\end{bmatrix}
$$

My test code was as follows:

```

int main()
{
  int n = 3;
  double A[9] = {67, 8, 4, 8, 90, 8, 4, 8, 83};
  double v[n];
  double x0[n];
  for (int i = 0; i < n; i++)
  {
    x0[i] = 1.0;
  }
  double mu = 77.0;
  for (int i = 0; i < n; i++)
    A[i * n + i] -= mu;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%f ", A[i * n + j]);
    }
    printf("\n");
  }
  printf("\n\n");
  double lambda = inversePowerMethod(n, A, v, x0, 1e-8, 100) + mu;
  printf("lambda = %f\n", lambda);
  printf("v = [");
  for (int i = 0; i < n; i++)
  {
    printf("%f", v[i]);
    if (i < n - 1)
      printf(", ");
  }
  printf("]\n");
  return 0;
}
```

I also added $\mu$ to the resulting eigenvalue. The result was:

```
lambda = 77.837045
v = [-0.079899, -0.617906, 1.019346]
```

This is the correct result, as the middle eigenvalue of $A$ is $77.837045$ and the corresponding eigenvector is $[-0.079899, -0.617906, 1.019346]$. This eigenvalue is different from the ones computed using the Inverse Power Method and the Power Method, which were $64.363547$ and $97.7995$, respectively.

## Task 3: Finding More Eigenvalues

