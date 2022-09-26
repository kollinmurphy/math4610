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
