from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np

import greedy3D


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#loc = tick.MultipleLocator(base = 0.4)
#ax.xaxis.set_major_locator(loc)
#ax.zaxis.set_major_locator(loc)
#ax.yaxis.set_major_locator(loc)

p = input("Protein: ")
times = int(input("How many times do you want to run the algorithm?"))

high_score = -1
scores = []
freqs = []
for i in range(0, times):
    seq, _, score = greedy3D.hill(len(p), p)
    print(score)
    if score not in scores:
        scores.append(score)
        freqs.append(1)
    else:
        freqs[scores.index(score)] += 1
    if score > high_score:
        high_score = score
        high_seq = seq

print("Score: -" + str(high_score))
z = []
y = []
x = []
for i in range(1, len(p) + 1):
    loc = np.where(high_seq == i)
    z.append(loc[0][0])
    y.append(loc[1][0])
    x.append(loc[2][0])
    ax.scatter(x[-1], y[-1], z[-1], s = 80, marker = r'$' + p[i - 1] + '$')

ax.set_xlim(min(x), max(x))
ax.set_ylim(min(y), max(y))
ax.set_zlim(min(z), max(z))

#ax.axis('off')

ax.plot_wireframe(x, y, z, color = 'yellow')

plt.show()
