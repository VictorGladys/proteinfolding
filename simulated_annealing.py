from includes.iterative_framework import *
import random

def anneal(max_iters, p, chance_threshold, T=None, scoref=None):
    if T == None:
        T = lambda i, maxi: i/maxi
    if scoref == None:
        scoref = init_score(p)
    seq, pos = init_grid(p, len(p))
    best_score = 0

    def continue(i, score):
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
                chance = random.randint(1, 100)
                if chance <= chance_threshold:
                    if this_score > best_score:
                        scores = [option]
                        best_score = this_score
                    elif this_score == best_score:
                        scores.append(option)
                else:
                    scores.append(option)
        seq, pos = random.choice(scores)
    return seq, pos, scoref(seq)

if __name__ == '__main__':
    iters = 10
    chance_threshhold = 80
    p = input("Protein: ")
    seq, pos, score = hill(iters, p, chance_threshhold)
    print(seq)
    print(score)
