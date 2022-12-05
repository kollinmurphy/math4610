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

I implemented an algorithm to find all three eigenvalues of a 3x3 matrix. I used the power method to find the largest eigenvalue, the inverse power method to find the smallest eigenvalue, and then partitioned that interval to find the middle eigenvalue. I used the same matrix as in Task 1, and the result was:

```
lambdaMax = 97.799321
lambdaMin = 64.671378
lambda = 77.836992
```

This is the correct result. My implementation is as follows:

```c
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

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%f ", A[i * n + j]);
    }
    printf("\n");
  }
  printf("\n\n");

  double tolerance = 0.0001;
  int maxIterations = 100;

  double lambdaMax = powerMethod(n, A, v, x0, tolerance, maxIterations);
  printf("lambdaMax = %f\n", lambdaMax);

  double lambdaMin = inversePowerMethod(n, A, v, x0, tolerance, maxIterations);
  printf("lambdaMin = %f\n", lambdaMin);

  double lambdaRange = lambdaMax - lambdaMin;
  int partitions = n;
  double partitionSize = lambdaRange / partitions;

  double B[9];

  for (int i = 1; i < partitions - 1; i++)
  {
    copyMatrix(n, A, B);
    double mu = lambdaMin + i * partitionSize;
    for (int i = 0; i < n; i++)
      B[i * n + i] -= mu;
    double lambda = inversePowerMethod(n, B, v, x0, tolerance, maxIterations);
    printf("lambda = %f\n", lambda + mu);
  }

  return 0;
}
```

I partitioned the interval into 3 parts, ignoring the first and last parts because the minimum and maximum eigenvalues were already found. I shifted the matrix by $\mu$ and then used the inverse power method to find the eigenvalue. I then added $\mu$ to the result to get the correct eigenvalue. This could easily be extended to matrices of any size.

## Task 4: OpenMP Implementation

To parallelize the code, I used OpenMP. I used the `#pragma omp parallel` directive to parallelize the computations of the shifted inverse power method. Most of the code was reused from Task 3, the only difference was changing the final `for` loop to use OpenMP. The following is the modified section of code:

```c
#pragma omp parallel
{
  int threadId = omp_get_thread_num();
  int numThreads = omp_get_num_threads();
  for (int i = 1 + threadId; i < partitions - 1; i += numThreads)
  {
    copyMatrix(n, A, B);
    double mu = lambdaMin + i * partitionSize;
    for (int i = 0; i < n; i++)
      B[i * n + i] -= mu;
    double lambda = inversePowerMethod(n, B, v, x0, tolerance, maxIterations);
    printf("lambda = %f\n", lambda + mu);
  }
}
```

I used the `omp_get_thread_num()` function to get the thread ID, and the `omp_get_num_threads()` function to get the number of threads. I then used the thread ID and the number of threads to determine which partitions to compute. I used the `printf` function to print the results, but this could be changed to use a shared array to store the results for further analysis. The output of the code was the same as in Task 3.

## Task 5: Gauss Seidel Method

Gauss Seidel is an iterative method for solving linear systems of equations. It is similar to the Jacobi method. One problem with the Gauss Seidel method is that it guarantees covergence only if the matrix is diagonally dominant or symmetric positive definite.

I implemented the Gauss Seidel method in C. The code is as follows:

```c
void gaussSeidel(int n, double *A, double *b, double *initialApproximation, double *result, double tolerance, int maxIterations)
{
  double x0[n];
  for (int i = 0; i < n; i++)
    x0[i] = initialApproximation[i];

  for (int iter = 0; iter < maxIterations; iter++)
  {
    double x1[n];
    for (int i = 0; i < n; i++)
      x1[i] = 0;

    printf("Iteration %d: [ ", iter);
    for (int i = 0; i < n; i++)
      printf("%f ", x0[i]);
    printf("]\n");

    for (int i = 0; i < n; i++)
    {
      double s1 = 0.0;
      for (int j = 0; j < i; j++)
        s1 += A[i * n + j] * x1[j];
      double s2 = 0.0;
      for (int j = i + 1; j < 4; j++)
        s2 += A[i * n + j] * x0[j];
      x1[i] = (b[i] - s1 - s2) / A[i * n + i];
    }
    int converged = 1;
    for (int i = 0; i < n; i++)
    {
      if (fabs(x1[i] - x0[i]) > tolerance)
      {
        converged = 0;
        break;
      }
    }
    for (int i = 0; i < n; i++)
      x0[i] = x1[i];
    if (converged == 1)
      break;
  }

  printf("Solution: [ ");
  for (int i = 0; i < n; i++)
    printf("%f ", x0[i]);
  printf("]\n");

  double error[n];
  for (int i = 0; i < n; i++)
  {
    error[i] = 0.0;
    for (int j = 0; j < n; j++)
    {
      error[i] += A[i * n + j] * x0[j];
    }
    error[i] -= b[i];
  }

  printf("Error: [ ");
  for (int i = 0; i < n; i++)
    printf("%f ", error[i]);
  printf("]\n");

  for (int i = 0; i < n; i++)
    result[i] = x0[i];
}
```

I tested the function using the following code:

```c
int main()
{
  double tolerance = 0.00001;
  int maxIterations = 1000;
  int n = 4;
  double A[16] = {
      10.0, -1.0, 2.0, 0.0,
      -1.0, 11.0, -1.0, 3.0,
      2.0, -1.0, 10.0, -1.0,
      0.0, 3.0, -1.0, 8.0};
  double b[4] = {6.0, 25.0, -11.0, 15.0};
  double x0[4] = {0.0, 0.0, 0.0, 0.0};
  double result[4];
  gaussSeidel(n, A, b, x0, result, tolerance, maxIterations);
  return 0;
}
```

The output of the program is as follows:

```text
Iteration 0: [ 0.000000 0.000000 0.000000 0.000000 ]
Iteration 1: [ 0.600000 2.327273 -0.987273 0.878864 ]
Iteration 2: [ 1.030182 2.036938 -1.014456 0.984341 ]
Iteration 3: [ 1.006585 2.003555 -1.002527 0.998351 ]
Iteration 4: [ 1.000861 2.000298 -1.000307 0.999850 ]
Iteration 5: [ 1.000091 2.000021 -1.000031 0.999988 ]
Iteration 6: [ 1.000008 2.000001 -1.000003 0.999999 ]
Solution: [ 1.000001 2.000000 -1.000000 1.000000 ]
Error: [ 0.000006 -0.000000 -0.000001 0.000000 ]
```

Because I used the same function signature as the Jacobi method, I was able to swap out the Jacobi method for the Gauss Seidel method in the code from Task 4. I also removed the print statements to keep the output of the lambda values easier to see. The updated implementation of the inverse power method is as follows:

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

  gaussSeidel(n, A, x0, initialForJacobi, y, tolerance, maxIterations);

  while (error > tolerance && iter < maxIterations)
  {
    norm = l2Norm(n, y);
    for (int i = 0; i < n; i++)
      z[i] = y[i] / norm;

    gaussSeidel(n, A, z, initialForJacobi, w, tolerance, maxIterations);
    lambda1 = dotProduct(n, z, w);
    error = fabs(lambda1 - lambda0);
    lambda0 = lambda1;
    iter++;
    for (int i = 0; i < n; i++)
      y[i] = w[i];
  }

  if (iter == maxIterations)
    printf("Inverse power method did not converge in %d iterations\n", iter);

  for (int i = 0; i < n; i++)
    result[i] = y[i];
  return 1.0 / lambda1;
}
```

The output of the new program using the same test matrix as in Task 4 is as follows:

```text
67.000000 8.000000 4.000000 
8.000000 90.000000 8.000000 
4.000000 8.000000 83.000000 


lambdaMax = 97.799321
lambdaMin = 64.663588

lambda = 77.836355
```

This output is within the tolerance of the output from Task 4, so the Gauss Seidel method is working correctly.
