import numpy as np

import simulated_annealing as sim
import includes.protein as protein

if __name__ == '__main__':
    prot = protein.Protein(None)
    n = prot.n

    # Get number of times to run algorithm
    times = int(input("How many times do you want to run the algorithm?"))

    # Initialize algorithm
    iters = 300
    high_score = -1
    scores = []
    freqs = []

    # Choice from various paths from 0 to 1 (y) over 0 tot iters (x)
    #f = sim.gen_exponentialT(iters, 0.01)
    f = sim.gen_linearT(iters)
    #f = sim.gen_oneT()
    #f = sim.gen_sigmoidT_mathv(iters)

    #Run algorithm n times
    for i in range(0, times):
        try:
            seq, _, score = sim.anneal(iters, prot.p, T=f)
        except KeyboardInterrupt:
            break

        print(score)
        print("Keer: ", i)

        # Log frequency of scores
        if score not in scores:
            scores.append(score)
            freqs.append(1)
        else:
            freqs[scores.index(score)] += 1

        # Log highest scores
        if score > high_score:
            high_score = score
            high_seq = seq
        

    # Visualize highest score
    np.set_printoptions(threshold=np.nan)
    print(high_seq)
    print("Score: -" + str(high_score))

    # Visualize frequencies
    print("Scores:     ", scores)
    print("Frequencies:", freqs)

    if input() == "s":
        for i in range(1, n + 1):
            print(np.where(high_seq != 0))

