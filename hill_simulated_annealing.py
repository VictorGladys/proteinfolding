from includes.iterative_framework import *
from includes.simulated_annealing import *
import random

# The approach we used with the greedy simulated annealing algorithm proved very
# slow, so we decided to cut down on the amount of times we'd run the score function.
# By doing this we implemented a true simulated annealing hillclimber.
# Every iteration we look at different options until we find one that scores high
# enough to be permitted by the temperature. We then take this option as the next
# starting point for our protein.
# This still means that at T is 1 we can only keep the current score or pick
# a higher score, whilst T = 0 permits any worse scoring option.
def anneal(max_iters, p, T=None, scoref=None):
    if T == None:
        T = gen_linearT(max_iters)
    if scoref == None:
        scoref = init_score(p)
    seq, pos = init_grid(p, len(p))
    global_best_score = 0
    global_best_fold = (seq, pos)
    all_perms = [(bendd, bendp) for bendd in [LEFT, RIGHT] for bendp in range(2, len(p))]
    fold = (seq, pos)

    for i in range(max_iters):
        score = 0
        thresh = T(i) * score

        # We search for options till we find one that is higher than the
        # threshold determined by the temperature. We hardcapped this
        # arbitrarily at 100 to prevent infinite loops.
        for _ in range(100):
            bendd, bendp = random.choice(all_perms)
            option = bend_part(bendd, bendp, seq, pos)
            if option[0] is None:
                continue

            score = scoref(option[0])
            if score >= thresh:
                seq, pos = option
                break
        else:
            print('Found only folds that make the protein worse.')

        if score > global_best_score:
            global_best_score = score
            global_best_fold = (seq, pos)
            print('new global best: ', global_best_score)

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
