import sys
import numpy as np
import time
from colorama import init, Back, Style
import dynamic_solver as dyn
init()

# This function shows the intermediary states of the dynamic programming solution
# in a somewhat more understandable way.
def grid_repr(q, startp, mid, endp):
    grid_repr = ''
    for i, line in enumerate(q):
        for j, char in enumerate(line):
            if  i == mid+1:
                grid_repr += Back.GREEN
            if (i, j) == (mid+1, endp):
                grid_repr += ' ' + Back.CYAN + str(char) + Style.RESET_ALL
            elif (i, j) == (startp, mid):
                grid_repr += ' ' + Back.YELLOW + str(char) + Style.RESET_ALL
            else:
                grid_repr += ' ' + str(char)
            if  i == mid+1:
                grid_repr += Style.RESET_ALL
        grid_repr += '\n'
    print(grid_repr)
    print('-'*50)
    time.sleep(0.5)

# This is the same function as the normal dynamic_solver algorithm, but passing
# an optional function argument to accomodate grid_repr in the normal algorithm
# would unnecessarily slow the algorithm down. Thusly we chose to accept not
# being able to refactor this code.
def solve(p):
    q    = [[0  for _ in range(len(p))] for _ in range(len(p))]
    fold = [[[] for _ in range(len(p))] for _ in range(len(p))]

    for mid in range(len(p) - 3, 0, -1): # == range(1, len(p) - 2) backwards
        for startp in range(0, mid):
            gains = []
            for endp in range(mid + 2, len(p)):
                val = dyn.Cprofit(p, startp, mid, endp)
                continues_from = q[mid+1][endp]

                gains.append(val + continues_from)

                print("Can gain by folding this way:", val,
                       "and continues from fold(s) that already summate: ", continues_from,
                       "\n which equals", val + continues_from)
                grid_repr(q, startp, mid, endp)
            idx = np.argmax(gains)
            q[startp][mid] = gains[idx]

            fold[startp][mid] += fold[mid][idx+mid+2] + [mid]
    grid_repr(q, startp, mid, endp)
    return q[0], fold[0]

# Initializes and runs the program
if __name__ == '__main__':
    p = input("Protein: ")

    q, fold = solve(p)
    maxval = np.argmax(q[1:len(p)]) + 1
    print("We read out the best value of the upper row of the finished grid, which is", maxval)

    print("To get this score, we have to fold after these indexes: ", fold[maxval])
