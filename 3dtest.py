from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')


ax1.axis('off')

x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
z = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

ax1.plot_wireframe(x, y, z, color = 'yellow')
ax1.scatter(x, y, z, s = 50, marker = r'$H$')

plt.show()
