from . import visualizer as vis
import numpy as np

# Protein object contains a list of aminoacids, that each have a type and location
class Protein(object):
    protein = []
    n = 0
    p = ""
    canv = None

    # Initialize
    def __init__(self, canv):
        self.p = input("Protein: ")
        self.n = len(self.p)
        self.canv = canv
        # Protein has to be at least 3 long to be relevant
        if self.n < 3:
            quit()

        # Protein has to be composed out of H and P and C only
        for i in self.p:
            if i != "H" and i != "P" and i != "C":
                quit()
            else:
                self.protein.append(self.Amino(i))

    def translateDynamic(self, folds):
        direction = 1
        last_w = self.n // 2
        last_h = self.n // 2
        skip_next = False

        for i in range(self.n):
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
            amino.drawAmino(self.protein, drawn, self.canv)
            drawn += 1
        vis.scoreFn(list(self.protein), 0)

    def translateIterative(self, seq):
        print(self.p)
        print(seq)
        last_w = self.n // 2
        last_h = np.ceil(self.n / 2)
        last_w_seq = self.n
        last_h_seq = self.n
        self.protein[0].assignPlace(last_w, last_h)

        for i in range(2, self.n + 1):
            if seq[last_h_seq][last_w_seq + 1] == i:
                last_w_seq += 1
                last_w += 1
            elif seq[last_h_seq][last_w_seq - 1] == i:
                last_w_seq -= 1
                last_w -= 1
            elif seq[last_h_seq + 1][last_w_seq] == i:
                last_h_seq += 1
                last_h += 1
            elif seq[last_h_seq - 1][last_w_seq] == i:
                last_h_seq -= 1
                last_h -= 1

            self.protein[i - 1].assignPlace(last_w, last_h)

        print("Score: -{}".format(vis.scoreFn(list(self.protein), 0)))

        drawn = 0
        for amino in self.protein:
            amino.drawAmino(self.protein, drawn, self.canv)
            drawn += 1
        vis.scoreFn(list(self.protein), 0)

    class Amino(object):
        prevAmino = None

        # Initialize
        def __init__(self, type):
            self.type = type

        def draw_line(self, canv, a, b, c, d):
            canv.create_line(vis.r(self.loc_w * 25 + 25 / a),
                               vis.r(self.loc_h * 25 + 25 / b),
                               vis.r(self.prevAmino.loc_w * 25 + 25 / c),
                               vis.r(self.prevAmino.loc_h * 25 + 25 / d), fill="red")

        # Draw an aminoacid
        def drawAmino(self, protein, drawn, canv):
            self.prevAmino = protein[drawn - 1]

            # Draw letter
            canv.create_text(vis.r(self.loc_w * 25 + 25 / 2), vis.r(self.loc_h * 25 + 25 / 2), text = self.type, font = vis.r(25), fill="red")


            # Connect with line to previous aminoacid
            if drawn > 0:
                if self.loc_h < self.prevAmino.loc_h:
                    self.draw_line(canv, 2, 1.5, 2, 3)
                elif self.loc_h > self.prevAmino.loc_h:
                    self.draw_line(canv, 2, 3, 2, 1.5)
                elif self.loc_w < self.prevAmino.loc_w:
                    self.draw_line(canv, 1.5, 2, 3, 2)
                else:
                    self.draw_line(canv, 3, 2, 1.5, 2)

        # Set location of aminoacid and request to draw
        def assignPlace(self, loc_w, loc_h):
            self.loc_w = loc_w
            self.loc_h = loc_h
