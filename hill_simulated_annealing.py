from includes.iterative_framework import *
import random

def gen_oneT():
    def oneT():
        return 1
    return oneT

def gen_linearT(maxi):
    def linearT(i):
        return i/maxi
    return linearT

def gen_exponentialT(maxi, b):
    g = (1/b)**(1/maxi)
    def exponentialT(i):
        return b*g**i
    return exponentialT

# f:         y0 + (y1 - y0)
#    -----------------------------
#                        x0+x1
#     (1 + exp(-v *(i -  ----- )
#                          2
#
# We try different approaches to calculate v, the
# average steepness of the slope
def gen_sigmoidT_mathv(maxi, eps=10e-4):
    x1 = maxi/2
    v = 2*np.log((1/eps) - 1) / maxi
    def sigmoidT(i):
        return 1 / (1 + np.exp(-v * (i - x1)))        
    return sigmoidT
                        
def gen_sigmoidT_tryv(maxi, eps=10e-6):
    x1 = maxi/2
    v = 1
    def sigmoidT(i):
        return 1 / (1 + np.exp(-v * (i - x1)))
    
    while not np.abs(    sigmoidT(0   )) < eps or\
          not np.abs(1 - sigmoidT(maxi)) < eps:
        v /= 2
        print(v)
        
    return sigmoidT
    

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
