import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("mesh.dat")
i = data[:, 0]
x = data[:, 1]

plt.plot(x, np.zeros_like(x), "o", markersize=3)
plt.xlabel("x")
plt.title("Mesh points (n = {})".format(len(x)))
plt.yticks([])
plt.grid(True, axis="x")
plt.show()
