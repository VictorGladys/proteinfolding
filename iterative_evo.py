from includes.iterative_framework import *
import random

# Big question: How do we deal with invalid bends? Two ways:
# 1) We ignore them, this would be akin to 'inactivated DNA'
# 2) We ignore those proteins. Risk of creating generations with only invalid
#    proteins.


# Two different ways of making offspring out of the previous generation:

# This method generates a new generation based on picking parents randomly from
# a skewed binomial curve
def generate_offspring1(idxs, gen_size, prev_gen, p_len):
        skewed = np.random.binomial(len(idxs), 0.9, gen_size//3)
        offspr = mix(prev_gen[idxs[skewed]], [random.randint(1, p_len-2)])
        return np.array(offspr)

# This method generates offspring at a fixed ratio, but generates a too large
# next generation
def generate_offspring2(idxs, gen_size, prev_gen, p_len):
        top = prev_gen[idxs[-top_n:]]
        more = (4*gen_size)//7
        less = gen_size - more
        select_offspr = mix(top, range(1, p_len-2))
        shared_offspr = mix(prev_gen, [p_len//2])

        return np.array(random.sample(select_offspr, more)+
                            random.sample(shared_offspr, less))

# Provides a method of creating children of every parent combination at every
# specified cutoffpoint.
def mix(proteins, cutoffpoint):
    new_proteins = []
    for i in cutoffpoint:
        for prota in proteins:
            for protb in proteins:
                new_protein = np.ndarray((len(p)-2,2))
                new_protein[:i] = prota[:i]
                new_protein[i:] = protb[i:]
                new_proteins.append(new_protein)
    return new_proteins

# Provides a random protein of a given length that is not necessarily a legal fold.
def random_protein(p_len):
    return [(random.choice([LEFT, RIGHT, NONE]), i) for i in range(2, p_len)]

def evo(gens, gen_size, top_n, p):
    seq, pos = init_grid(p, len(p))
    score = init_score(p)
    best_score = 0

    this_gen = np.array([random_protein(len(p)) for _ in range(gen_size)])

    for i in range(gens):
        prev_gen = this_gen
        # First we find the best proteins in the previous generation...
        scores = []
        for parent in prev_gen:
            # should construct a protein from a list of [bendd, bendp] pairs
            new_seq, _ = bend_all(parent, seq, pos)
            if new_seq is None:
                scores.append(-1)
            else:
                scores.append(score(new_seq))

        rand = np.random.random(len(scores))

        idxs = np.lexsort((rand, scores))

        # ...and then we base the new generation on them.
        # we need a mechanism to imporve diversity/uniqueness
        # otherwise we get stuck on local minima too easily
        this_gen = generate_offspring1(idxs, gen_size, prev_gen, len(p))
        print(i)


    new_seq, _ = bend_all(this_gen[-1], seq, pos)
    print(new_seq)
    print(score(new_seq))

if __name__ == '__main__':

    gens = 100
    gen_size = 100

    top_n = 5
    p = input("Protein: ")

    evo(gens, gen_size, top_n, p)
