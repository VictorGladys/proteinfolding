import sys
import numpy as np
import time
from colorama import init, Back, Style
init()

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

if __name__ == '__main__':
    p = input("Protein: ")
    q = []
    fold = [[] for _ in range(len(p))]
    for mid in range(1, len(p)):
        q.append([])
        for n in range(len(p)):
            q[mid-1].append(0)

    for mid in range(len(p) - 3, 0, -1):
        for startp in range(0, mid):
            gains = [0]
            for endp in range(mid + 2, len(p)):
                val = Cprofit(p, startp, mid, endp)
                continues_from = q[mid+1][endp]

                gains.append(val + continues_from)

                print("Can gain by folding this way:", val,
                       "and continues from fold(s) that already summate: ", continues_from,
                       "\n which equals", val + continues_from)
                grid_repr(q, startp, mid, endp)
            q[startp][mid] = max(gains)
maxval = max(q[0][1:len(p)])
print("Readout best value of the upper row of the finished grid:\n which is", maxval)
grid_repr(q, startp, mid, endp)

