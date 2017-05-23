from includes.iterative_framework import *
import random

def linearT(i, maxi):
    return i/maxi

def anneal(max_iters, p, T=None, scoref=None):
    if T == None:
        T = linearT
    if scoref == None:
        scoref = init_score(p)
    seq, pos = init_grid(p, len(p))
    best_score = 0
    global_best_score = 0

    for i in range(max_iters):
        scores = [[best_score, [seq, pos]]]
        for bendd in [LEFT, RIGHT]:
            # We ignore 'bends' in the first and last positions
            # because they are nonsensical
            # additionaly, counting starts from one, because of the way our array
            # is implemented
            for bendp in range(2, len(p)):
                option = bend_part(bendd, bendp, seq, pos)
                if option[0] is None:
                    continue

                this_score = scoref(option[0])
                print(option)
                scores.append([this_score, list(option)])
                print(option)
                return
        print(scores[0])
        scores = np.array(scores[0])
        thresh = T(i, max_iters) * best_score
        coords = np.where(scores[:, 0] >= thresh)
        print(coords)
        best_score, (seq, pos) = random.choice(scores[coords])

        print("HOI1")
        tmp = np.argmax(scores[:, 0])
        if scores[tmp, 0] > global_best_score:
            global_best_score = scores[tmp, 0]
            global_best_fold = scores[tmp, 1]
    
    return global_best_fold[0], global_best_fold[1], global_best_score

if __name__ == '__main__':
    iters = 10
    p = input("Protein: ")
    seq, pos, score = anneal(iters, p)
    print(seq)
    print(score)
