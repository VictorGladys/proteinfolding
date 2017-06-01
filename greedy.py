from includes.iterative_framework import *
import random

# Despite being called hill, this functions implements a greedy search through
# 2D state space. For every iteration, it calculates the score and its legality,
# of any additional fold to the protein. Then, with every score known,
# it accepts one of the best scoring folds at random.
def hill(iters, p, scoref=None):
    seq, pos = init_grid(p, len(p))
    if scoref == None:
        scoref = init_score(p)
    best_score = 0

    for i in range(iters):
        scores = [(seq, pos)]
        for bendd in [LEFT, RIGHT]:
            # We ignore 'bends' in the first and last positions
            # because they are nonsensical
            # additionaly, counting starts from one, because of the way our array
            # is implemented
            for bendp in range(2, len(p)):
                option = bend_part(bendd, bendp, seq, pos)
                if option[0] is None:
                    continue

                # If an option comes along that scores equally as good as the
                # current best scoring option, it will be remembered as well.
                # If we see a better scoring option, all previous options will
                # be disregarded.
                this_score = scoref(option[0])
                if this_score > best_score:
                    scores = [option]
                    best_score = this_score
                elif this_score == best_score:
                    scores.append(option)
        # Next we choose randomly from the options with the best score.
        seq, pos = random.choice(scores)
    return seq, pos, scoref(seq)

if __name__ == '__main__':
    iters = 100
    p = input("Protein: ")
    seq, pos, score = hill(iters, p)
    print(seq)
    print(score)
