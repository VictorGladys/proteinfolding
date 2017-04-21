import numpy as np

LEFT  = 1
RIGHT = 2

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
def bend_part(bend_dir, pos, seq_grid, pos_grid):
    start_coord = np.where(seq_grid == pos)
    coords = [[], []]
    new_coords = [[], []]

    # find coordinates of following aminoacids
    while 1:
        pos += 1
        coord = np.where(seq_grid == pos)
        if coord[0].size == 0:
            break
        else:
            coords[0].append(coord[0])
            coords[1].append(coord[1])

    # translate to origin
    new_coords[0] = np.array(coords[0]) - start_coord[0] # y values
    new_coords[1] = np.array(coords[1]) - start_coord[1] # x values

    # turn coordinates
    direction = 1 if bend_dir == RIGHT else -1
    temp = new_coords[0]
    new_coords[0] = direction*new_coords[1]
    new_coords[1] = temp

    # retranslate back to cutoffpoint
    new_coords[0] = np.array(new_coords[0]) + start_coord[0] # y values
    new_coords[1] = np.array(new_coords[1]) + start_coord[1] # x values


    # check for incorrect fold
    pos_grid[coords] = 0
    pos_grid[new_coords] += 1

    if np.where(pos_grid > 1)[0].size > 0:
        print("incorrect fold occured")
        return 0

    seq_grid[new_coords] = seq_grid[coords]
    seq_grid[coords] = 0

    print(seq_grid)
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
    p = input("Protein: ")

    seq_grid, pos_grid = init_grid(p, len(p))
    bend_part(RIGHT, 2, seq_grid, pos_grid)
    bend_part(RIGHT, 4, seq_grid, pos_grid)
    bend_part(RIGHT, 6, seq_grid, pos_grid)

