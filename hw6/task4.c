#include <stdio.h>
#include <math.h>

double factorial(int n)
{
  if (n == 0)
  {
    return 1;
  }
  return n * factorial(n - 1);
}

int main()
{
  double sum = 0;
  int n = 1000;
  for (int i = n; i >= 0; i--)
  {
    sum += 1 / factorial(i);
  }
  printf("e = %.16f\n", sum);
  return 0;
}
