from tkinter import *
import matplotlib.pyplot as plt
from numpy import array

import simulated_annealing as sim
import includes.protein as protein
import includes.visualizer as vis

if __name__ == '__main__':
    # Initialize window
    root=Tk()
    root.attributes("-topmost", True)
    root.title("Proteam Power")

    # Globals

    canvas = Canvas(root)

    prot = protein.Protein(canvas)
    n = prot.n
    times = int(input("How many times do you want to run the algorithm?"))

    # Initialize canvas
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

    iters = 1000
    high_score = -1
    scores = []
    freqs = []
    #f = sim.gen_exponentialT(iters, 0.01)
    #f = sim.gen_linearT(iters)
    #f = sim.gen_oneT()
    f = sim.gen_sigmoidT_mathv(iters)
    for i in range(0, times):
        seq, _, score = sim.anneal(iters, prot.p, T=f)
        print(score)
        if score not in scores:
            scores.append(score)
            freqs.append(1)
        else:
            freqs[scores.index(score)] += 1
        if score > high_score:
            high_score = score
            high_seq = seq
            
    prot.translateIterativeHill(high_seq)

    npscores = array(scores)
    npfreqs = array(freqs)
    plt.scatter(npscores, npfreqs)
    plt.show()
    

    root.mainloop()
