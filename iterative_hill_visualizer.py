from tkinter import *

import iterative_hill as itr
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

    q, fold = dyn.solve(prot.p)
    maxval = np.argmax(q[1:n]) + 1
    prot.translateDynamic(fold[maxval])

    root.mainloop()


