## Vesa Frijling
## visualisatieProtein
## 04-11-2017
##
## Dit programma visualiseert hoe hydrofobe en hydrofile aminozuren samen en proteine vormen

from tkinter import *

# Adjusts to different screen sizes, change r to change canvas size
def r(number):
    ratio = 0.8
    if int(number):
        return int(number * ratio)
    else:
        return -1

def scoreFn(to_check, s):
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
        return scoreFn(to_check, s)
    else:
        return s

# Aminoacid object contains a type and location    
class Amino (object):
    
    # Initialize
    def __init__(self, type):
        self.type = type

    # Draw an aminoacid    
    def drawAmino(self):
        global canvas, drawn, protein

        # Draw letter
        id = canvas.create_text(r(self.loc_w * grid_size + grid_size / 2), r(self.loc_h * grid_size + grid_size / 2), text = self.type, font = r(grid_size), fill="red")
        letters.append(id)

        # Connect with line to previous aminoacid
        if drawn > 0:
            if self.loc_h < protein[drawn - 1].loc_h:
                id = canvas.create_line(r(self.loc_w * grid_size + grid_size / 2), r(self.loc_h * grid_size + grid_size / 1.5),
                                   r(protein[drawn - 1].loc_w * grid_size + grid_size / 2), r(protein[drawn - 1].loc_h * grid_size + grid_size / 3), fill="red")
            elif self.loc_h > protein[drawn - 1].loc_h:
                id = canvas.create_line(r(self.loc_w * grid_size + grid_size / 2), r(self.loc_h * grid_size + grid_size / 3),
                                   r(protein[drawn - 1].loc_w * grid_size + grid_size / 2), r(protein[drawn - 1].loc_h * grid_size + grid_size / 1.5), fill="red")
            elif self.loc_w < protein[drawn - 1].loc_w:
                id = canvas.create_line(r(self.loc_w * grid_size + grid_size / 1.5), r(self.loc_h * grid_size + grid_size / 2),
                                   r(protein[drawn - 1].loc_w * grid_size + grid_size / 3), r(protein[drawn - 1].loc_h * grid_size + grid_size / 2), fill="red")
            else:
                id = canvas.create_line(r(self.loc_w * grid_size + grid_size / 3), r(self.loc_h * grid_size + grid_size / 2),
                                   r(protein[drawn - 1].loc_w * grid_size + grid_size / 1.5), r(protein[drawn - 1].loc_h * grid_size + grid_size / 2), fill="red")
        lines.append(id)        

        # Adjust amount of aminoacids drawn
        drawn += 1

        if drawn == n:
            print(scoreFn(list(protein), 0))

    # Set location of aminoacid and request to draw
    def assignPlace(self, loc_w, loc_h):
        self.loc_w = loc_w
        self.loc_h = loc_h
        

# Validate whether clicked location is valid for next aminoacid in protein
def validatePlace(event, loc_w, loc_h):
    # Check if any aminacids left to process
    if drawn == n:
        return

    # Check if another aminoacid already holds intended location
    for i in range(drawn - 1):
        if protein[i].loc_h == loc_h and protein[i].loc_w == loc_w:
            return

    # Check if in spot directly next to previous aminoacid
    if abs(loc_w - protein[drawn - 1].loc_w) ==  1 and loc_h == protein[drawn - 1].loc_h:
        protein[drawn].assignPlace(loc_w, loc_h)
        protein[drawn].drawAmino()
    elif abs(loc_h - protein[drawn - 1].loc_h) ==  1 and loc_w == protein[drawn - 1].loc_w:
        protein[drawn].assignPlace(loc_w, loc_h)
        protein[drawn].drawAmino()

# Delete previous aminoacid except for first two        
def delete():
    global drawn
    if drawn > 2:
        canvas.delete(lines[-1])
        del lines[-1]
        
        canvas.delete(letters[-1])
        del letters[-1]
        
        drawn -= 1

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
        if i != "H" and i != "P" and i != "C":
            quit()
        else:
            protein.append(Amino(i))

# Helper for on click event for grid item
# When cell clicked calls validatePlace
def tag_helper(loc_w, loc_h):
        return lambda event:validatePlace(event, loc_w, loc_h)

# Setup grid
def createGrid():
    for i in range(2 * n + 1):
        for j in range(2 * n + 1):
            id = canvas.create_rectangle(r(j * grid_size), r(i * grid_size), r((j + 1) * grid_size), r((i + 1) * grid_size), fill='white')
            canvas.tag_bind(id, "<Button-1>", tag_helper(j, i))
    canvas.pack()
    b = Button(root, text="Back", command=delete)
    b.pack(side=BOTTOM)

if __name__ == '__main__':
    # Initialize window
    root=Tk()
    root.attributes("-topmost", True)
    root.title("Proteam Power")

    # Globals
    grid_size = 30
    drawn = 0
    lines = []
    letters = []
    protein = []
    
    setupProtein()

    # Initialize canvas
    canvas = Canvas(root, width = r((2 * n + 1) * grid_size), height = r((2 * n + 1) * grid_size))
    createGrid()

    # Place first and second aminoacid
    protein[drawn].assignPlace(n, n)
    protein[drawn].drawAmino()
    protein[drawn].assignPlace(n, n - 1)
    protein[drawn].drawAmino()
    root.mainloop()
