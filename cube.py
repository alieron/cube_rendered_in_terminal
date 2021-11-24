from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import matplotlib.pyplot as plt

from math import sqrt

fig = plt.figure()
ax = plt.axes(projection="3d")


# data = np.array(
#     [[1, 0, -1 / sqrt(2)],
#      [-1, 0, -1 / sqrt(2)],
#      [0, 1, 1 / sqrt(2)],
#      [0, -1, 1 / sqrt(2)]])

data = []

for i in np.arange(-1, 1.4, 0.4):
    for j in np.arange(-1, 1.4, 0.4):
        for k in np.arange(-1, 1.4, 0.4):
            if i % 1 == 0 or j % 1 == 0 or k % 1 == 0:
                data.append([i, j, k])

# print(data)

data = np.array(data)

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

ax.scatter3D(x, y, z, color='red')

# connections = [[0, 1, 2, 3, 0], [0, 2], [1, 3]]

# for digit in connections:
#     X = [x[a] for a in digit]
#     Y = [y[a] for a in digit]
#     Z = [z[a] for a in digit]

#     ax.plot3D(X, Y, Z, color="green")

ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))

plt.show()
