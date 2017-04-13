import sys
import numpy as np

# We simply count the amount of bridges formed by folding behind mid,
# only counting from and to start and end
# mid    ---    mid+1
#  |             |
# mid-it      mid+1+it
def profit(protein, start, mid, end):
    shortest = min(mid - start, end - (mid+1))
    bridges = 0
    for it in range(1, shortest+1):
        if protein[mid - it] == 'H' and protein[mid + 1 + it] == 'H':
            bridges += 1
    return bridges


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
                # is it really [mid ---- 1]?
                val = profit(p, startp, mid, endp) + q[mid+1][endp]
                gains.append(val)
            idx = np.argmax(gains)
            q[startp][mid] = gains[idx]
        # op welke manier moet je deze indices nu opslaan?
        fold[idx].append(mid)
q = q[0]
maxval = np.argmax(q[1:len(p)]) + 1
print("Folds after indexes: ", fold[maxval])
print(fold)
print(q)
print("Score: -", q[maxval])
