import numpy as np
import copy

E = ( 1, 0)
W = (-1, 0)
N = ( 0, 1)
S = ( 0,-1)
ds = set([E, W, N, S])

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

    bridges = 0
    if grid[pos[1]][pos[0]] == 'H':
        for neighd in neighbordirs:
            try:
                neighpos = pos + neighd
                if grid[neighpos[1]][neighpos[0]] == 'H':
                    bridges += 1
            except:
                pass
    if grid[pos[1]][pos[0]] == 'C':
        for neighd in neighbordirs:
            try:
                neighpos = pos + neighd
                if grid[neighpos[1]][neighpos[0]] == 'C':
                    bridges += 5
            except:
                pass

    if rest[1:] == '':
        #print('\n'.join(' '.join(' ' if char == '' else char for char in line) for line in grid))
        #print('-'*50)
        return bridges

    return bridges + max([0] + [next(newd, pos, copy.deepcopy(grid), rest[1:])
                for newd in neighbordirs
                    if grid[pos[1] + newd[1]][pos[0] + newd[0]] == ''])




if __name__ == '__main__':
    p = input("Protein: ")
    length = len(p)
    pos = np.array([max(0, length-4), length-1])
    width = length + max(pos[0], 0)
    height = 2 * (length - 1) + 1

    grid = [['' for _ in range(width)] for _ in range(height)]
    grid[pos[1]][pos[0]] = p[0]

    print(next(E, pos, grid, p[1:]))
