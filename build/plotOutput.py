
import numpy as np
from matplotlib import pyplot as plt
import time

E = np.genfromtxt('output.csv', delimiter=',')

print("E.shape : ", E.shape)
Nt = E.shape[0]

plt.ion()
plt.figure()
for i in range(Nt):
  plt.plot(E[i,:], 'b')
  plt.pause(0.05)
  plt.clf()
  
plt.show()


