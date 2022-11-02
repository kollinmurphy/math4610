def explicitEulerLogistic(alpha, beta, p0, n=100, h=1.0):
  def f(p):
    return alpha * p - beta * (p * p)

  t0 = 0
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
