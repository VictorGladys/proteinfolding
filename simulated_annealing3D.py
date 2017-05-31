from includes.iterative_framework import *
import random

def gen_oneT():
    def oneT(i):
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
        scoref = init_score(p, is_3d=True)
    seq, pos = init_grid(p, len(p))
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

        scores = np.array(scores)
        thresh = T(i) * best_score
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
