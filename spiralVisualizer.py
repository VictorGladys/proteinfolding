## Vesa Frijling
## visualisatieProtein
## 04-11-2017
##
## Dit programma visualiseert hoe hydrofobe en hydrofile aminozuren samen en proteine vormen in een spiraalvorm

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
            if j != to_check[1] and (j.type == "H"  or j.type == "C") \
                 and ((check.loc_h == j.loc_h and abs(check.loc_w - j.loc_w) == 1) \
                 or (check.loc_w == j.loc_w and abs(check.loc_h - j.loc_h) == 1)):
                s += 1
    if check.type == "C" and len(to_check) > 1:
        for j in to_check:
            if j != to_check[1] and j.type == "H" \
                 and ((check.loc_h == j.loc_h and abs(check.loc_w - j.loc_w) == 1) \
                 or (check.loc_w == j.loc_w and abs(check.loc_h - j.loc_h) == 1)):
                s += 1
            if j != to_check[1] and j.type == "C" \
                 and ((check.loc_h == j.loc_h and abs(check.loc_w - j.loc_w) == 1) \
                 or (check.loc_w == j.loc_w and abs(check.loc_h - j.loc_h) == 1)):
                s += 5
                         
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

        
def spiral(up, result):
    x = 0
    result[0].assignPlace(n // 2, n // 2)
    counter = 1
    last_h = n // 2
    last_w = n // 2
    
    while counter < n:
        for i in range(up + x):
            if counter < n:
                result[counter].assignPlace(last_w, last_h - 1)
                last_h -= 1
                counter += 1
                
        x += 1
        
        for i in range(x):
            if counter < n:
                result[counter].assignPlace(last_w + 1, last_h)
                last_w += 1
                counter += 1
                
        for i in range(up + x):
            if counter < n:
                result[counter].assignPlace(last_w, last_h + 1)
                last_h += 1
                counter += 1
                                
        x += 1

        for i in range(x):
            if counter < n:
                result[counter].assignPlace(last_w - 1, last_h)
                last_w -= 1
                counter += 1
                
    scoreFn(list(result), 0)
    return [score, list(result)]

# Get protein
def setupProtein():
    global n
    
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

# Setup grid
def createGrid():
    for i in range(n + 1):
        for j in range(n + 1):
            canvas.create_rectangle(r(j * grid_size), r(i * grid_size), r((j + 1) * grid_size), r((i + 1) * grid_size), fill='white')
            
    canvas.pack()

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

    for i in range(round((n - 4) // 2)):
        results.append(spiral(i + 1, list(protein)))
        #copy = list(protein)
        #copy.reverse()
        #results.append(spiral(i + 1, copy))

    top = 0
    index = 0
    for i in range(len(results)):
        if results[i][0] > top:
            top = results[i][0]
            index = i
    print(top)
    for i in spiral(index + 1, list(protein))[1]:
        i.drawAmino()
    
    root.mainloop()
