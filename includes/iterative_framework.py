#
# This functions as a framework for our iterative algorithms. Include this with:
# import iterative_framework
#
# Use it like demo'ed in this codes main

import numpy as np
import copy

LEFT  = 0
NONE  = 1
RIGHT = 2
UP    = 3
DOWN  = 4

def cost(aminoA, aminoB):
    charsum = ord(aminoA) + ord(aminoB)
    if charsum == 144 or charsum == 139:    # H + H or H + C
        return 1
    elif charsum == 134:                    # C + C
        return 5
    return 0

# seq_grid:     protein:
#  00000          HPPHHH
#  06120
#  05430
#  00000
#
# We check the neighbors of i=1 if they contain a pos higher than i+1. One of
# the neighbors are 4 and 6, so we check if they add to the score. They do, so
# our score becomes 2. As we increment i, after this initial scoring we don't
# encounter anymore neighbor_pos_val > this_pos_val+1, so the score stays 2.
def init_score(p, is_3d=False):
    def neigh2d(seq_grid, coord):
        N = seq_grid[coord[0]-1, coord[1]  ][0]-1
        S = seq_grid[coord[0]+1, coord[1]  ][0]-1
        O = seq_grid[coord[0]  , coord[1]+1][0]-1
        W = seq_grid[coord[0]  , coord[1]-1][0]-1
        return [N, S, O, W]

    def neigh3d(seq_grid, coord):
        N = seq_grid[coord[0]  , coord[1]-1, coord[2]  ][0]-1
        S = seq_grid[coord[0]  , coord[1]+1, coord[2]  ][0]-1
        O = seq_grid[coord[0]  , coord[1]  , coord[2]+1][0]-1
        W = seq_grid[coord[0]  , coord[1]  , coord[2]-1][0]-1
        U = seq_grid[coord[0]-1, coord[1]  , coord[2]  ][0]-1
        D = seq_grid[coord[0]+1, coord[1]  , coord[2]  ][0]-1
        return [N, S, O, W, U, D]

    neighfun = neigh3d if is_3d else neigh2d

    def score(seq_grid):
        s = 0
        # -1 does not lose information because we only care if a pos is bigger
        # also: this way we never go out of bounds with the directions
        for pos in range(len(p)-1):
            coord = np.where(seq_grid == pos+1)

            neighs = neighfun(seq_grid, coord)

            # pos+1 because we don't want to count sequential aminos
            for neigh in neighs:
                if neigh > pos+1:
                    s += cost(p[pos], p[neigh])
        return s
    return score

#                       N
#           1        (0, -2)       -1
#                    (0, -1)
# W  (-2, 0) (-1, 0) (0,  0) (1, 0)  (2, 0)  E
#                    (0,  1)
#          -1        (0,  2)        1
#                       S
# The way we rotate coordinates is thusly:
# We can see a bend a a point in our protein as a linear rotation. Since
# a rotation in linear algebra can only rotate over the origin, we translate
# the coordinates to the origin. Our point x now sits in the above picture at
# point 0, 0.
# Next we may notice that moving from W to S is simply exhanging te values,
# but moving from S to W is exhanging the values but also multiplying by -1,
# I've added the multiplicator s needed to correctly rotate to above picture in
# the corners.
#
# We may decribe the correct value of s as so:
# Turning right is positive and left negative, except it is switched if we move
# along the x-axis prior to turning/ if we do not move along the y-axis.
def rotate_coords(coords, start_coord, prev_coord, bend_dir):
    new_coords = [[], []]

    # translate to origin
    new_coords[0] = np.array(coords[0]) - start_coord[0] # y values
    new_coords[1] = np.array(coords[1]) - start_coord[1] # x values

    # turn coordinates
    s = 1 if bend_dir == RIGHT else -1

    ## next we check if we move along the y axis by taking the difference
    ## between this and the previous y-coordinate. 0 means no, 1 or -1 yes
    s *= -1 if start_coord[0] - prev_coord[0] else 1
    temp = new_coords[0]
    new_coords[0] = new_coords[1] * s
    new_coords[1] = temp          * s

    # retranslate back to cutoffpoint
    new_coords[0] = np.array(new_coords[0]) + start_coord[0] # y values
    new_coords[1] = np.array(new_coords[1]) + start_coord[1] # x values

    return new_coords

#                              U
#                          (-2, 0, 0)
#                          (-1, 0, 0)
#
#                                N
#                  1         (0, 0, -2)       -1
#                           (0, 0, -1)
# W  (0, -2, 0) (0, -1, 0) (0, 0,  0) (0, 1, 0)  (0, 2, 0)  E
#                         (0, 0,  1)
#              -1        (0, 0,  2)        1
#                            S
#
#                          (1, 0, 0)
#                          (2, 0, 0)
#                              D
# N->U switch 2 to 0
# E->U switch 1 to 0
# S->U switch 2 to 0  *-1
# W->U switch 1 to 0  *-1
#
# N->D switch 2 to 0  *-1
# E->D switch 1 to 0  *-1
# S->D switch 2 to 0
# W->D switch 1 to 0
def rotate_coords3d(coords, start_coord, prev_coord, bend_dir):
    new_coords = [[], [], []]

    # translate to origin
    new_coords[0] = np.array(coords[0]) - start_coord[0] # z values
    new_coords[1] = np.array(coords[1]) - start_coord[1] # y values
    new_coords[2] = np.array(coords[2]) - start_coord[2] # x values

    # turn coordinates
    if bend_dir in [UP, DOWN]:
        s = 1 if bend_dir == DOWN else -1
        ## next we check the directionality along the x and y axis
        ## with this we can find out the correct value for s, and wether we
        ## need to switch the z and x, or z and y coordinates
        y_directionality = start_coord[1] - prev_coord[1]
        x_directionality = start_coord[2] - prev_coord[2]

        c = 1 if y_directionality  else 2
        s *= y_directionality if y_directionality else x_directionality

        temp = new_coords[0]
        new_coords[0] = new_coords[c] * s
        new_coords[c] = temp          * s
    elif bend_dir in [RIGHT, LEFT]:
        s = 1 if bend_dir == RIGHT else -1
        s *= -1 if start_coord[0] - prev_coord[0] else 1
        temp = new_coords[1]
        new_coords[1] = new_coords[2] * s
        new_coords[2] = temp          * s

    # retranslate back to cutoffpoint
    new_coords[0] = np.array(new_coords[0]) + start_coord[0] # z values
    new_coords[1] = np.array(new_coords[1]) + start_coord[1] # y values
    new_coords[2] = np.array(new_coords[2]) + start_coord[2] # x values

    return new_coords

def find_relevant_coords(seq_grid, bend_pos):
    coords = [[], []]

    # find coordinates of following aminoacids
    while 1:
        bend_pos += 1
        coord = np.where(seq_grid == bend_pos)
        if coord[0].size == 0:
            break
        else:
            coords[0].append(coord[0])
            coords[1].append(coord[1])
    return coords

def find_relevant_coords3d(seq_grid, bend_pos):
    coords = [[], [], []]

    # find coordinates of following aminoacids
    while 1:
        bend_pos += 1
        coord = np.where(seq_grid == bend_pos)
        # check if the bend_pos exists/ has coordinates
        if coord[0].size == 0:
            break
        else:
            coords[0].append(coord[0])
            coords[1].append(coord[1])
            coords[2].append(coord[2])
    return coords

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
def bend_part(bend_dir, bend_pos, orig_seq_grid, orig_pos_grid, is_3d=False):
    seq_grid = copy.deepcopy(orig_seq_grid)
    pos_grid = copy.deepcopy(orig_pos_grid)

    prev_coord = np.where(seq_grid == bend_pos-1)
    start_coord = np.where(seq_grid == bend_pos)

    if not is_3d:
        coords = find_relevant_coords(seq_grid, bend_pos)
        new_coords = rotate_coords(coords, start_coord, prev_coord, bend_dir)
    else:
        coords = find_relevant_coords3d(seq_grid, bend_pos)
        new_coords = rotate_coords3d(coords, start_coord, prev_coord, bend_dir)

    # check for incorrect fold
    pos_grid[coords] = 0
    pos_grid[new_coords] += 1

    if np.where(pos_grid > 1)[0].size > 0:
        #print("incorrect fold occured, skipped")
        return None, None

    tmp_vals = seq_grid[coords]
    seq_grid[coords] = 0
    seq_grid[new_coords] = tmp_vals

    return seq_grid, pos_grid


def bend_all(bends, orig_seq, orig_pos):
    seq = copy.deepcopy(orig_seq)
    pos = copy.deepcopy(orig_pos)

    for bend_dir, bend_pos in bends:
        if bend_dir == NONE:
            continue
        seq, pos = bend_part(bend_dir, bend_pos, seq, pos)
        if seq is None:
            return None, None

    return seq, pos


def init_grid(p, l):
    seq_grid = np.ndarray((2*l, 2*l), dtype=int)
    seq_grid[:, :] = 0
    pos_grid = np.ndarray((2*l, 2*l), dtype=int)
    pos_grid[:, :] = 0

    seq_grid[(l, l + np.arange(l))] = np.arange(l)+1
    pos_grid[(l, l + np.arange(l))] = 1

    return seq_grid, pos_grid

def init_space(p, l):
    seq_space = np.ndarray((2*l, 2*l, 2*l), dtype=int)
    seq_space[:, :, :] = 0
    pos_space = np.ndarray((2*l, 2*l, 2*l), dtype=int)
    pos_space[:, :, :] = 0

    seq_space[(l, l + np.arange(l))] = np.arange(l)+1
    pos_space[(l, l + np.arange(l))] = 1

    return seq_grid, pos_grid


if __name__ == '__main__':
    p = "HPPHHPHPPH"
    seq, pos = init_grid(p, len(p))
    score = init_score(p)

    seq, pos = bend_part(LEFT, 2, seq, pos)
    print(seq)
    print("Score: ", score(seq))
    seq, pos = bend_part(LEFT, 4, seq, pos)
    print(seq)
    print("Score: ", score(seq))
    seq, pos = bend_part(LEFT, 6, seq, pos)
    print(seq)
    print("Score: ", score(seq))
    seq, pos = bend_part(LEFT, 9, seq, pos)
    print(seq)
    print("Score: ", score(seq))
