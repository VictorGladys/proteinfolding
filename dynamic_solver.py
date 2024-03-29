import sys
import numpy as np

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

# This algorithm calculates every possible standalone fold- and while doing so
# adds up the best scoring compatible standalone fold that came prior to the
# current standalone fold. By doing this we limit ourself to a accordion fold,
# limiting our maximum score, but also reduce the complexity from O(3^n-2)
# to O(n^4). (n^3 because of three for loops, times the cost function which runs in O(n))
def solve(p):
    q    = [[0  for _ in range(len(p))] for _ in range(len(p))]
    fold = [[[] for _ in range(len(p))] for _ in range(len(p))]

    for mid in range(len(p) - 3, 0, -1): # == range(1, len(p) - 2) backwards
        for startp in range(0, mid):
            gains = []
            for endp in range(mid + 2, len(p)):
                val = Cprofit(p, startp, mid, endp)
                continues_from = q[mid+1][endp]

                gains.append(val + continues_from)

            # We use argmax to be able both work with the score as well as the
            # storing the thusly scoring protein in a way that we can find it
            # after we know the highest score
            idx = np.argmax(gains)
            q[startp][mid] = gains[idx]

            fold[startp][mid] += fold[mid][idx+mid+2] + [mid]
    return q[0], fold[0]

# Initializes and runs the program
if __name__ == '__main__':
    p = input("Protein: ")

    q, fold = solve(p)
    maxval = np.argmax(q[1:len(p)]) + 1

    print("Folds after indexes: ", fold[maxval])
    print("Score: -{}".format(q[maxval]))
