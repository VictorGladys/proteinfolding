from includes.iterative_framework import *
from includes.simulated_annealing import *
import random

# The approach we used with the greedy algorithm continues into our simulated
# annealing algorithm- we calculate the score of every fold that is possible in
# one step from the current protein, then choose randomly out of the subset of
# options that scores higher than the threshold described by the temperature times
# the score of the current protein. This means that at T is 1 we can only
# keep the current score or pick a higher score, whilst T = 0 permits any worse
# scoring option.
def anneal(max_iters, p, T=None, scoref=None):
    if T == None:
        T = gen_linearT(max_iters)
    if scoref == None:
        scoref = init_score(p, is_3d=True)
    seq, pos = init_grid3d(p, len(p))
    best_score = 0
    global_best_score = 0
    global_best_fold = None

    i = 0
    while i < max_iters:
        scores = [best_score]
        folds = [(seq, pos)]
        for bendd in [LEFT, RIGHT, UP, DOWN]:
            # We ignore 'bends' in the first and last positions
            # because they are nonsensical
            # additionaly, counting starts from one, because of the way our array
            # is implemented
            for bendp in range(2, len(p)):
                option = bend_part(bendd, bendp, seq, pos, is_3d=True)
                if option[0] is None:
                    continue

                this_score = scoref(option[0])
                scores.append(this_score)
                folds.append(option)

        # We calculate the threshold based on a temperature function that we pass
        # to the anneal function. T has been initialised so that it only needs an
        # 0 <= i <= max_iters to find the correct temperature.
        scores = np.array(scores)
        thresh = T(i) * best_score

        # Here we then take a random coordinate of one of the scores that is permitted
        # by the temperature, and then use it to find the corresponding folded protein.
        coords = np.where(scores >= thresh)
        coord = random.choice(coords[0])

        best_score = scores[coord]
        (seq, pos) = folds[coord]

        tmp = np.argmax(scores)
        if scores[tmp] > global_best_score:
            global_best_score = scores[tmp]
            global_best_fold = folds[tmp]
            print('new global best: ', global_best_score, ', at iteration: ', i)
            i = int(i*0.75)

        i += 1

    return global_best_fold[0], global_best_fold[1], global_best_score

if __name__ == '__main__':
    iters = 100
    p = input("Protein: ")
    seq, pos, score = anneal(iters, p)
    print(seq)
    print(score)

    f = gen_exponentialT(iters, 0.01)
    seq, pos, score = anneal(iters, p, T=f)
    print(seq)
    print(score)
