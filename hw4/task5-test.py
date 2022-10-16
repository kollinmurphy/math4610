import matplotlib.pyplot as plt
from hw4.task5 import explicitEulerLogistic

tvals, pvals = explicitEulerLogistic(2.0, 0.0005, 10)
print(tvals, pvals)
plt.title('Test 3')
plt.xlabel('t')
plt.ylabel('P(t)')
plt.plot(tvals, pvals)
plt.show()
