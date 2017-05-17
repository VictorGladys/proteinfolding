from includes.iterative_framework import *
import random

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

                this_score = scoref(option[0])
                if this_score > best_score:
                    scores = [option]
                    best_score = this_score
                elif this_score == best_score:
                    scores.append(option)
        seq, pos = random.choice(scores)
    return seq, pos, scoref(seq)

if __name__ == '__main__':
    iters = 100
    p = input("Protein: ")
    seq, pos, score = hill(iters, p)
    print(seq)
    print(score)