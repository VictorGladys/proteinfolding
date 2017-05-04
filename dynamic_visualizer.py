import sys
import numpy as np
from tkinter import *

import includes/dynamic_solver as dyn
import includes/visualizer as vis

# Protein object contains a list of aminoacids, that each have a type and location
class Protein (object):
    protein = []
    n = 0
    p = ""

    # Initialize
    def __init__(self):
        self.p = input("Protein: ")
        self.n = len(self.p)

        # Protein has to be at least 3 long to be relevant
        if self.n < 3:
            quit()

        # Protein has to be composed out of H and P and C only
        for i in self.p:
            if i != "H" and i != "P" and i != "C":
                quit()
            else:
                self.protein.append(self.Amino(i))

    def translate(self, folds):
        direction = 1
        last_w = n // 2
        last_h = n // 2
        skip_next = False

        for i in range(n):
            if skip_next == False:
                self.protein[i].assignPlace(last_w, last_h + direction)
                last_h += direction
                if i in folds:
                    self.protein[i + 1].assignPlace(last_w + 1, last_h)
                    skip_next = True
                    last_w += 1
                    direction *= -1
            else:
                skip_next = False

        print("Score: -{}".format(vis.scoreFn(list(self.protein), 0)))

        drawn = 0
        for amino in self.protein:
            amino.drawAmino(self.protein, drawn)
            drawn += 1
        vis.scoreFn(list(self.protein), 0)

    class Amino(object):
        prevAmino = None

        # Initialize
        def __init__(self, type):
            self.type = type

        def draw_line(self, canvas, a, b, c, d):
            canvas.create_line(vis.r(self.loc_w * grid_size + grid_size / a),
                               vis.r(self.loc_h * grid_size + grid_size / b),
                               vis.r(self.prevAmino.loc_w * grid_size + grid_size / c),
                               vis.r(self.prevAmino.loc_h * grid_size + grid_size / d), fill="red")

        # Draw an aminoacid
        def drawAmino(self, protein, drawn):
            global canvas
            self.prevAmino = protein[drawn - 1]

            # Draw letter
            canvas.create_text(vis.r(self.loc_w * grid_size + grid_size / 2), vis.r(self.loc_h * grid_size + grid_size / 2), text = self.type, font = vis.r(grid_size), fill="red")


            # Connect with line to previous aminoacid
            if drawn > 0:
                if self.loc_h < self.prevAmino.loc_h:
                    self.draw_line(canvas, 2, 1.5, 2, 3)
                elif self.loc_h > self.prevAmino.loc_h:
                    self.draw_line(canvas, 2, 3, 2, 1.5)
                elif self.loc_w < self.prevAmino.loc_w:
                    self.draw_line(canvas, 1.5, 2, 3, 2)
                else:
                    self.draw_line(canvas, 3, 2, 1.5, 2)

        # Set location of aminoacid and request to draw
        def assignPlace(self, loc_w, loc_h):
            self.loc_w = loc_w
            self.loc_h = loc_h

# Setup grid
def createGrid():
    for i in range(n + 1):
        for j in range(n + 1):
            canvas.create_rectangle(vis.r(j * grid_size), vis.r(i * grid_size),
                    vis.r((j + 1) * grid_size), vis.r((i + 1) * grid_size), fill='white')
    canvas.pack()

if __name__ == '__main__':
    # Initialize window
    root=Tk()
    root.attributes("-topmost", True)
    root.title("Proteam Power")

    # Globals
    grid_size = 25

    prot = Protein()
    n = prot.n

    # Initialize canvas
    canvas = Canvas(root, width = vis.r((n + 1) * grid_size),
                          height = vis.r((n + 1) * grid_size),
                          scrollregion=(0, 0, vis.r((n + 1) * grid_size),
                          vis.r((n + 1) * grid_size)))

    hbar = Scrollbar(root,orient=HORIZONTAL)
    hbar.pack(side=BOTTOM,fill=X)
    hbar.config(command=canvas.xview)

    vbar = Scrollbar(root,orient=VERTICAL)
    vbar.pack(side=RIGHT,fill=Y)
    vbar.config(command=canvas.yview)

    canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
    canvas.pack(side=LEFT,expand=True,fill=BOTH)

    createGrid()

    q, fold = dyn.solve(prot.p)
    maxval = np.argmax(q[1:n]) + 1
    prot.translate(fold[maxval])

    root.mainloop()
