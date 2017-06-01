from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np

import hill_simulated_annealing3D as hlsim
from includes.simulated_annealing import *

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

p = input("Protein: ")
times = int(input("How many times do you want to run the algorithm?"))
iters = 3000

# Choice from various paths from 0 to 1 (y) over 0 tot iters (x)
#f = gen_exponentialT(iters, 0.01)
#f = gen_linearT(iters)
#f = gen_oneT()
f = gen_sigmoidT_mathv(iters)

high_score = -1
scores = []
freqs = []
for i in range(0, times):
    print("Time: ", i)
    try:
        seq, _, score = hlsim.anneal(iters, p, T=f)
    except KeyboardInterrupt:
        break

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
