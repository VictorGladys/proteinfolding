#
# This functions as a framework for our iterative algorithms. Include this with:
# import iterative_framework
#
# Use it like demo'ed in this codes main

import numpy as np
import copy

LEFT  = 1
RIGHT = 2

def cost(aminoA, aminoB):
    charsum = ord(aminoA) + ord(aminoB)
    if charsum == 144 or charsum == 139:    # H + H or H + C
        return 1
    elif charsum == 134:                    # C + C
        return 5
    return 0

def score(p, seq_grid):
    s = 0
    # -1 does not lose information because we only care if a pos is bigger
    # also: this way we never go out of bounds with the directions
    for pos in range(len(p)-1):
        coord = np.where(seq_grid == pos+1)

        N = seq_grid[coord[0]+1, coord[1]  ][0]-1
        S = seq_grid[coord[0]-1, coord[1]  ][0]-1
        O = seq_grid[coord[0]  , coord[1]+1][0]-1
        W = seq_grid[coord[0]  , coord[1]-1][0]-1

        # pos+1 because we don't want to count sequential aminos
        if N > pos+1:
            s += cost(p[pos], p[N])
        elif S > pos+1:
            s += cost(p[pos], p[S])
        elif O > pos+1:
            s += cost(p[pos], p[O])
        elif W > pos+1:
            s += cost(p[pos], p[W])
    return s


def rotate_coords(coords, start_coord, prev_coord, bend_dir):
    new_coords = [[], []]

    # translate to origin
    new_coords[0] = np.array(coords[0]) - start_coord[0] # y values
    new_coords[1] = np.array(coords[1]) - start_coord[1] # x values

    # turn coordinates
    s = 1 if bend_dir == RIGHT else -1
    s *= -1 if start_coord[0] - prev_coord[0] else 1
    temp = new_coords[0]
    new_coords[0] = new_coords[1] * s
    new_coords[1] = temp          * s

    # retranslate back to cutoffpoint
    new_coords[0] = np.array(new_coords[0]) + start_coord[0] # y values
    new_coords[1] = np.array(new_coords[1]) + start_coord[1] # x values

    return new_coords

# We want to bend after pos, so:
# 1 - 2 - 3
# and
# pos is 2, bend_dir is DOWN
#  |
#  v
#
# 1 - 2
#     |
#     3
def bend_part(bend_dir, pos, orig_seq_grid, orig_pos_grid):
    seq_grid = copy.deepcopy(orig_seq_grid)
    pos_grid = copy.deepcopy(orig_pos_grid)

    prev_coord = np.where(seq_grid == pos-1)
    start_coord = np.where(seq_grid == pos)
    coords = [[], []]

    # find coordinates of following aminoacids
    while 1:
        pos += 1
        coord = np.where(seq_grid == pos)
        if coord[0].size == 0:
            break
        else:
            coords[0].append(coord[0])
            coords[1].append(coord[1])

    new_coords = rotate_coords(coords, start_coord, prev_coord, bend_dir)

    # check for incorrect fold
    pos_grid[coords] = 0
    pos_grid[new_coords] += 1

    if np.where(pos_grid > 1)[0].size > 0:
        print("incorrect fold occured, skipped")
        return orig_seq_grid, orig_pos_grid

    seq_grid[new_coords] = seq_grid[coords]
    seq_grid[coords] = 0

    return seq_grid, pos_grid

def init_grid(p, l):
    seq_grid = np.ndarray((2*l, 2*l), dtype=int)
    seq_grid[:, :] = 0
    pos_grid = np.ndarray((2*l, 2*l), dtype=int)
    pos_grid[:, :] = 0

    seq_grid[(l, l + np.arange(l))] = np.arange(l)+1
    pos_grid[(l, l + np.arange(l))] = 1

    return seq_grid, pos_grid


if __name__ == '__main__':
    p = "HPPHHPHPPH"

    seq, pos = init_grid(p, len(p))
    seq, pos = bend_part(LEFT, 2, seq, pos)
    print(seq)
    print("Score: ", score(p, seq))
    seq, pos = bend_part(LEFT, 4, seq, pos)
    print(seq)
    print("Score: ", score(p, seq))
    seq, pos = bend_part(LEFT, 6, seq, pos)
    print(seq)
    print("Score: ", score(p, seq))
    seq, pos = bend_part(LEFT, 9, seq, pos)
    print(seq)
    print("Score: ", score(p, seq))

