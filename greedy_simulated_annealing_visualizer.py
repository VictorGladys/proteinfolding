from tkinter import *
import matplotlib.pyplot as plt
from numpy import array

import greedy_simulated_annealing as sim
import includes.protein as protein
import includes.visualizer as vis
from includes.simulated_annealing import *

if __name__ == '__main__':
    # Initialize window
    root=Tk()
    root.attributes("-topmost", True)
    root.title("Proteam Power")

    # Globals
    canvas = Canvas(root)

    prot = protein.Protein(canvas)
    n = prot.n

    # Get number of times to run algorithm
    times = int(input("How many times do you want to run the algorithm?"))

    # Initialize canvas (scrollbar)
    canvas.config(width = vis.r((n + 1) * 25),
                          height = vis.r((n + 1) * 25),
                          scrollregion=(0, 0, vis.r((n + 1) * 25),
                          vis.r((n + 1) * 25)))

    hbar = Scrollbar(root,orient=HORIZONTAL)
    hbar.pack(side=BOTTOM,fill=X)
    hbar.config(command=canvas.xview)

    vbar = Scrollbar(root,orient=VERTICAL)
    vbar.pack(side=RIGHT,fill=Y)
    vbar.config(command=canvas.yview)

    canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
    canvas.pack(side=LEFT,expand=True,fill=BOTH)

    vis.createGrid(n, canvas)

    # Initialize algorithm
    iters = 3000
    high_score = -1
    scores = []
    freqs = []

    # Choice from various paths from 0 to 1 (y) over 0 tot iters (x)
    #f = gen_exponentialT(iters, 0.01)
    f = gen_linearT(iters)
    #f = gen_oneT()
    #f = gen_sigmoidT_mathv(iters)

    #Run algorithm n times
    for i in range(0, times):
        try:
            seq, _, score = sim.anneal(iters, prot.p, T=f)
        except KeyboardInterrupt:
            break

        print(score)
        print("Keer: ", i)

        # Log frequency of scores
        if score not in scores:
            scores.append(score)
            freqs.append(1)
        else:
            freqs[scores.index(score)] += 1

        # Log highest scores
        if score > high_score:
            high_score = score
            high_seq = seq

    # Visualize highest score
    prot.translateIterativeHill(high_seq)

    # Visualize frequencies
    npscores = array(scores)
    npfreqs = array(freqs)
    plt.scatter(npscores, npfreqs)
    plt.show()

    root.mainloop()
