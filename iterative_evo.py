from includes/iterative_framework import *
import random

# Big question: How do we deal with invalid bends? Two ways:
# 1) We ignore them, this would be akin to 'inactivated DNA'
# 2) We ignore those proteins. Risk of creating generations with only invalid
#    proteins.

def random_lengths(gen_size, p):
    return (np.ceil(np.clip((np.random.randn(gen_size)+1.2) *p/3, 0, p-1)))

def evo(gens, gen_size, p):
    seq, pos = init_grid(p, len(p))
    score = init_score(p)
    best_score = 0
    all_perms = [(bendd, bendp) for bendd in [LEFT, RIGHT] for bendp in range(2, len(p))]

# make generation_size amount of children with max(len(p)) amount of bends

    this_gen = [random.sample(all_perms, l) for l in random_lengths(gen_size, p)]
    print(this_gen)

    for _ in range(gens):
        prev_gen = this_gen
        # First we find the best proteins in the previous generation...
        scores = []
        for parent in prev_gen:
            # should construct a protein from a list of [bendd, bendp] pairs
            seq, pos = bend_all(parent)
            if seq is None:
                scores.append(-1)
            else:
                scores.append(score(seq))
        # we have to find the best 10% or something (this will lead to an unpredictable amount of children)
        # or we take the best 5 (even if 5th place is tied) this will lead to the same amount of children every time
        scores = np.array(scores)
        idxs = np.where(scores > max(scores)*0.9)

        # ...and then we base the new generation on them.
        this_gen = "mix"(prev_gen[idxs])
        # "mix"() will have to do the following: take each of the previous gen proteins
        # 1) take another one and take halve of each parents bends
        # 2) by chance replace or add one bend with a random new bend
        # 3) randomly remove a bend
        #
        # This means for x picked parents, we will get x! new children if we have each protein pair up with each other protein
        # or x children if each protein only pairs up with 1 other one

    print(seq)
    print(score(seq))

if __name__ is '__main__':
    gens = 100
    gen_size = 10
    p = input("Protein: ")

    evo(gens, gen_size, p)
