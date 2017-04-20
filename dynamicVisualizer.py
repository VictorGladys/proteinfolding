import sys
import numpy as np
from tkinter import *

# Adjusts to different screen sizes, change r to change canvas size
def r(number):
    ratio = 1
    if int(number):
        return int(number * ratio)
    else:
        return -1

def scoreFn(to_check, s):
    global score
    check = to_check[0]
    if check.type == "H" and len(to_check) > 1:
        for j in to_check:
            if j != to_check[1] and j.type == "H" \
                 and ((check.loc_h == j.loc_h and abs(check.loc_w - j.loc_w) == 1) \
                 or (check.loc_w == j.loc_w and abs(check.loc_h - j.loc_h) == 1)):
                s += 1
    del to_check[to_check.index(check)]
    if to_check != []:
        scoreFn(to_check, s)
    else:
        score = s

# Aminoacid object contains a type and location    
class Amino (object):
    
    # Initialize
    def __init__(self, type):
        self.type = type

    # Draw an aminoacid    
    def drawAmino(self):
        global canvas, drawn, protein

        # Draw letter
        canvas.create_text(r(self.loc_w * grid_size + grid_size / 2), r(self.loc_h * grid_size + grid_size / 2), text = self.type, font = r(grid_size), fill="red")

        # Connect with line to previous aminoacid
        if drawn > 0:
            if self.loc_h < protein[drawn - 1].loc_h:
                canvas.create_line(r(self.loc_w * grid_size + grid_size / 2), r(self.loc_h * grid_size + grid_size / 1.5),
                                   r(protein[drawn - 1].loc_w * grid_size + grid_size / 2), r(protein[drawn - 1].loc_h * grid_size + grid_size / 3), fill="red")
            elif self.loc_h > protein[drawn - 1].loc_h:
                canvas.create_line(r(self.loc_w * grid_size + grid_size / 2), r(self.loc_h * grid_size + grid_size / 3),
                                   r(protein[drawn - 1].loc_w * grid_size + grid_size / 2), r(protein[drawn - 1].loc_h * grid_size + grid_size / 1.5), fill="red")
            elif self.loc_w < protein[drawn - 1].loc_w:
                canvas.create_line(r(self.loc_w * grid_size + grid_size / 1.5), r(self.loc_h * grid_size + grid_size / 2),
                                   r(protein[drawn - 1].loc_w * grid_size + grid_size / 3), r(protein[drawn - 1].loc_h * grid_size + grid_size / 2), fill="red")
            else:
                canvas.create_line(r(self.loc_w * grid_size + grid_size / 3), r(self.loc_h * grid_size + grid_size / 2),
                                   r(protein[drawn - 1].loc_w * grid_size + grid_size / 1.5), r(protein[drawn - 1].loc_h * grid_size + grid_size / 2), fill="red")       

        # Adjust amount of aminoacids drawn
        drawn += 1

        if drawn == n:
            scoreFn(list(protein), 0)

    # Set location of aminoacid and request to draw
    def assignPlace(self, loc_w, loc_h):
        self.loc_w = loc_w
        self.loc_h = loc_h

# Setup grid
def createGrid():
    for i in range(n + 1):
        for j in range(n + 1):
            canvas.create_rectangle(r(j * grid_size), r(i * grid_size), r((j + 1) * grid_size), r((i + 1) * grid_size), fill='white')
            
    canvas.pack()

# We simply count the amount of bridges formed by folding behind mid,
# only counting from and to start and end
# mid    ---    mid+1
#  |             |
# mid-it      mid+1+it
def Hprofit(protein, start, mid, end):
    shortest = min(mid - start, end - (mid+1))
    bridges = 0
    for it in range(1, shortest+1):
        if protein[mid - it] == 'H' and protein[mid + 1 + it] == 'H':
            bridges += 1
    return bridges

def Cprofit(protein, start, mid, end):
    shortest = min(mid - start, end - (mid+1))
    bridges = 0
    for it in range(1, shortest+1):
        charsum = ord(protein[mid - it]) + ord(protein[mid + 1 + it])
        if charsum == 144:      # H + H
            bridges += 1
        elif charsum == 139:    # H + C
            bridges += 1
        elif charsum == 134:    # C + C
            bridges += 5
    return bridges

# Get protein
def setupProtein():
    global n, p
    
    p = input("Protein: ")
    n = len(p)

    # Protein has to be at least 3 long to be relevant
    if n < 3:
        quit()

    # Protein has to be composed out of H and P only
    for i in p:
        if i != "H" and i != "P":
            quit()
        else:
            protein.append(Amino(i))

def translate(folds):
    global protein
    
    direction = 1
    last_w = n // 2
    last_h = n // 2
    skip_next = False
    for i in range(n):
        if skip_next == False:
            protein[i].assignPlace(last_w, last_h + direction)
            last_h += direction
            if i in folds:
                protein[i + 1].assignPlace(last_w + 1, last_h)
                skip_next = True
                last_w += 1
                direction *= -1
        else:
            skip_next = False

    scoreFn(list(protein), 0)
    print("Score: -{}".format(score))
    
    for i in protein:
        i.drawAmino()

if __name__ == '__main__':
    # Initialize window
    root=Tk()
    root.attributes("-topmost", True)
    root.title("Proteam Power")
    
    # Globals
    grid_size = 25
    drawn = 0
    protein = []
    results = []
    p = ""
    
    setupProtein()

    # Initialize canvas
    canvas = Canvas(root, width = r((n + 1) * grid_size), height = r((n + 1) * grid_size), scrollregion=(0,0,r((n + 1) * grid_size),r((n + 1) * grid_size)))
    hbar=Scrollbar(root,orient=HORIZONTAL)
    hbar.pack(side=BOTTOM,fill=X)
    hbar.config(command=canvas.xview)
    vbar=Scrollbar(root,orient=VERTICAL)
    vbar.pack(side=RIGHT,fill=Y)
    vbar.config(command=canvas.yview)
    canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
    canvas.pack(side=LEFT,expand=True,fill=BOTH)
    createGrid()
    
    q    = [[0  for _ in range(len(p))] for _ in range(len(p))]
    fold = [[[] for _ in range(len(p))] for _ in range(len(p))]

    for mid in range(len(p) - 3, 0, -1): # == range(1, len(p) - 2) backwards
        for startp in range(0, mid):
            gains = [0]
            for endp in range(mid + 2, len(p)):
                val = Cprofit(p, startp, mid, endp)
                continues_from = q[mid+1][endp]

                gains.append(val + continues_from)

            idx = np.argmax(gains)
            q[startp][mid] = gains[idx]

            fold[startp][mid] += fold[mid][idx+1+mid] + [mid]
            #print('\n'.join('\t'.join(''.join(str(c) for c in char) for char in line) for line in fold))
            #print('-'*50)
    q = q[0]
    maxval = np.argmax(q[1:len(p)]) + 1
    translate(fold[0][maxval])
    
    
