import numpy as np
import copy

E = ( 1, 0)
W = (-1, 0)
N = ( 0, 1)
S = ( 0,-1)
ds = set([E, W, N, S])

# Recursively does a depth first search
def next(d, pos, grid, rest):
    possible_neighbordirs = (ds - set([(-d[0], -d[1])]))
    neighbordirs = set(possible_neighbordirs)
    for newd in possible_neighbordirs:
        newpos = pos + newd
        if pos[0] >= width or pos[1] >= height or pos[0] < 0 or pos[1] < 0:
            neighbordirs.remove(newd)
    d = np.array(d)
    pos = pos + d

    grid[pos[1]][pos[0]] = rest[0]

    # These two statements go through all neighbors and calculate the intermediary
    # score that the unfinished protein holds
    bridges = 0
    if grid[pos[1]][pos[0]] == 'H':
        for neighd in neighbordirs:
            try:
                neighpos = pos + neighd
                if grid[neighpos[1]][neighpos[0]] in ['H', 'C']:
                    bridges += 1
            except:
                pass
    elif grid[pos[1]][pos[0]] == 'C':
        for neighd in neighbordirs:
            try:
                neighpos = pos + neighd
                if grid[neighpos[1]][neighpos[0]] == 'H':
                    bridges += 1
                elif grid[neighpos[1]][neighpos[0]] == 'C':
                    bridges += 5
            except:
                pass

    # If we are at the outermost leaves of the searchtree, we return the curent
    # score
    if rest[1:] == '':
        return bridges

    # We return the highest score of all children spawned from this position
    return bridges + max([0] + [next(newd, pos, copy.deepcopy(grid), rest[1:])
                for newd in neighbordirs
                    if grid[pos[1] + newd[1]][pos[0] + newd[0]] == ''])

# Initializes and runs the program
if __name__ == '__main__':
    p = input("Protein: ")
    length = len(p)
    pos = np.array([max(0, length-4), length-1])
    width = length + max(pos[0], 0)
    height = 2 * (length - 1) + 1

    grid = np.array([['' for _ in range(width)] for _ in range(height)])
    grid[pos[1]][pos[0]] = p[0]

    print(next(E, pos, grid, p[1:]))
    input()
